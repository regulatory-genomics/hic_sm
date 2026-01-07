# Copied from https://github.com/mirnylab/distiller-sm/blob/master/_distiller_common.py
from pathlib import Path
import os
import sys
import yaml
import re
import numpy as np
import pandas as pd
import shlex
from peppy import Project
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("8.20.1")

module_name = "hic_pipeline"



# --- Config & PEP Loading ---

# Load PEP project to get sample metadata
pep_config_value = config.get("pep_config", os.path.join("pep", "project_config.yaml"))
# Resolve PEP path; if relative, treat it as relative to the repo root (one level above workflow/)
if os.path.isabs(pep_config_value):
    pep_config_path = pep_config_value
else:
    pep_config_path = os.path.abspath(os.path.join(workflow.basedir, "..", pep_config_value))

if not os.path.exists(pep_config_path):
    raise OSError(f"Project config file path does not exist: {pep_config_path}")

# Load PEP project
pep_project = Project(pep_config_path)

# Merge PEP path variables into config for easier access
if "paths" in pep_project.config:
    if "paths" not in config:
        config["paths"] = {}
    config["paths"].update(pep_project.config["paths"])

# Helper function to resolve path placeholders in config
def resolve_config_paths(cfg, paths_dict):
    """Recursively resolve path placeholders in config dictionary"""
    if isinstance(cfg, dict):
        return {k: resolve_config_paths(v, paths_dict) for k, v in cfg.items()}
    elif isinstance(cfg, list):
        return [resolve_config_paths(item, paths_dict) for item in cfg]
    elif isinstance(cfg, str) and "{" in cfg and "}" in cfg:
        for key, value in paths_dict.items():
            cfg = cfg.replace(f"{{{key}}}", value)
        return cfg
    else:
        return cfg

# Resolve path placeholders in config
# First, resolve PEP paths (for placeholders like {process_dir}, {database_dir})
paths_dict = {}
if "paths" in config:
    paths_dict.update(config["paths"])
    config = resolve_config_paths(config, paths_dict)

# Then resolve output_base_dir itself if it contains placeholders
if "output_base_dir" in config and isinstance(config["output_base_dir"], str) and "{" in config["output_base_dir"]:
    config["output_base_dir"] = resolve_config_paths(config["output_base_dir"], paths_dict)

# Finally, add output_base_dir to paths_dict and resolve {output_base_dir} placeholders in output paths
if "output_base_dir" in config:
    paths_dict["output_base_dir"] = config["output_base_dir"]
    config = resolve_config_paths(config, paths_dict)

# Get sample table from PEP
annot = pep_project.sample_table.copy()

# Explode list-valued columns if any
list_cols = [col for col in annot.columns if annot[col].apply(lambda x: isinstance(x, list)).any()]
if list_cols:
    # Explode all list columns simultaneously to preserve index alignment
    # This prevents Cartesian product issues when multiple columns contain lists
    # (e.g., when run=[1,2] and R1=[fileA,fileB] should stay paired as run=1->fileA, run=2->fileB)
    annot = annot.explode(list_cols, ignore_index=True)

# Remove duplicate (sample_name, run) pairs
if 'sample_name' in annot.columns and 'run' in annot.columns:
    annot = annot.drop_duplicates(subset=['sample_name', 'run'], keep='first').reset_index(drop=True)

# Ensure required columns exist
required_columns = {"sample_name", "run"}
missing_required = required_columns.difference(annot.columns)
if missing_required:
    raise ValueError(f"Sample table missing required columns: {', '.join(sorted(missing_required))}.")

# Ensure run column is integer for schema validation
annot['run'] = pd.to_numeric(annot['run'], errors='raise').astype(int)

# Ensure passqc column is integer if it exists
if 'passqc' in annot.columns:
    annot['passqc'] = pd.to_numeric(annot['passqc'], errors='coerce').fillna(0).astype(int)

# Normalize skip_ligation to boolean if present
if 'skip_ligation' in annot.columns:
    annot['skip_ligation'] = (
        annot['skip_ligation']
        .astype(str)
        .str.strip()
        .str.lower()
        .map({'true': True, 'false': False})
    )
    # Default to False for any unrecognized/NA values
    annot['skip_ligation'] = annot['skip_ligation'].fillna(False).astype(bool)

# Create a unique identifier for each run: sample_name_run (store run as string for downstream paths)
annot['sample_name'] = annot['sample_name'].astype(str)
annot['sample_run'] = annot['sample_name'] + '_run' + annot['run'].astype(str)

# Validate sample table
schema_path = os.path.join(workflow.basedir, "..", "schemas", "sample.schemas.yml")
if os.path.exists(schema_path):
    validate(annot, schema=schema_path)

# Convert run back to string for downstream path building
annot['run'] = annot['run'].astype(str)

# Check for duplicate sample_run identifiers
duplicates = annot[annot.duplicated(subset=['sample_run'], keep=False)]
if not duplicates.empty:
    sys.stderr.write("\nError: Duplicate (sample_name + run) pairs found in sample sheet!\n")
    sys.stderr.write("The following rows are not unique:\n")
    sys.stderr.write(str(duplicates[['sample_name', 'run', 'sample_run']]) + "\n\n")
    sys.exit(1)

# Set sample_run as index
annot = annot.set_index('sample_run')

# Get unique sample names (for downstream analysis)
samples = annot.reset_index().drop_duplicates(subset='sample_name', keep='first').set_index("sample_name").to_dict(orient="index")

# --- Path Definitions ---

# Set up output directories - define directly using output_base_dir
result_path = config.get("output_base_dir", "results")
outdir = result_path

# Genome configuration
assembly = config["input"]["genome"]["assembly_name"]
genome_path = config["input"]["genome"]["bwa_index_wildcard_path"].rstrip("*")
bowtie_index_path = config["input"]["genome"]["bowtie_index_path"].rstrip("*")
chromsizes_path = config["input"]["genome"]["chrom_sizes_path"]

# Define index targets based on mapper
MAPPER = config["map"]["mapper"]
if MAPPER == "bwa-mem":
    idx = multiext(
        genome_path,
        ".amb",
        ".ann",
        ".bwt",
        ".pac",
        ".sa",
    )
elif MAPPER == "bwa-mem2":
    idx = multiext(
        genome_path,
        ".0123",
        ".amb",
        ".ann",
        ".bwt.2bit.64",
        ".pac",
    )
elif MAPPER == "bowtie2":
    idx = multiext(
        bowtie_index_path,
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )
else:  # bwa-meme
    idx = multiext(
        genome_path,
        ".0123",
        ".amb",
        ".ann",
        ".pac",
        ".pos_packed",
        ".suffixarray_uint64",
        ".suffixarray_uint64_L0_PARAMETERS",
        ".suffixarray_uint64_L1_PARAMETERS",
        ".suffixarray_uint64_L2_PARAMETERS",
    )

# --- Sample Aggregation ---

# Build SAMPLE_FASTQS from annot - aggregate all runs per sample
SAMPLE_FASTQS = {}  # sample_name -> {"R1": [list], "R2": [list]}
SAMPLE_METADATA = {}  # sample_name -> metadata

for sample_run_id, row in annot.iterrows():
    sample_name = row['sample_name']
    
    # Get fastq paths (support both R1/R2 and fastq1/fastq2)
    if 'R1' in row and 'R2' in row:
        fastq1, fastq2 = row['R1'], row['R2']
    elif 'fastq1' in row and 'fastq2' in row:
        fastq1, fastq2 = row['fastq1'], row['fastq2']
    else:
        raise ValueError(f"Sample {sample_run_id} missing R1/R2 or fastq1/fastq2 columns")
    
    if sample_name not in SAMPLE_FASTQS:
        SAMPLE_FASTQS[sample_name] = {"R1": [], "R2": []}
        SAMPLE_METADATA[sample_name] = {}
    
    # Collect all FASTQ files for this sample
    SAMPLE_FASTQS[sample_name]["R1"].append(str(fastq1))
    SAMPLE_FASTQS[sample_name]["R2"].append(str(fastq2))
    
    # Store metadata (use first run's metadata or merge if needed)
    if not SAMPLE_METADATA[sample_name]:
        if 'ligation_site' in row and pd.notna(row['ligation_site']):
            SAMPLE_METADATA[sample_name]['ligation_site'] = str(row['ligation_site'])
        if 'skip_ligation' in row and pd.notna(row['skip_ligation']):
            SAMPLE_METADATA[sample_name]['skip_ligation'] = bool(row['skip_ligation'])

# Keep LIBRARY_RUN_FASTQS for backward compatibility with library_group logic
# (library = sample_name in current implementation)
LIBRARY_RUN_FASTQS = {sample: {"1": SAMPLE_FASTQS[sample]["R1"]} for sample in SAMPLE_FASTQS.keys()}

# Automatically create library_groups if all_group is True
if config["input"].get("all_group", False):
    all_samples = list(SAMPLE_FASTQS.keys())
    if "library_groups" not in config["input"]:
        config["input"]["library_groups"] = {}
    config["input"]["library_groups"]["all"] = all_samples
    print(f"Auto-generated 'all' group with {len(all_samples)} libraries: {all_samples}")

# --- Wildcard Constraints ---

# Get all sample_run identifiers for wildcard constraints
all_sample_runs = list(annot.index)

# Global constraints on the wildcards
wildcard_constraints:
    sample=f"({'|'.join([re.escape(s) for s in SAMPLE_FASTQS.keys()])})",
    library=f"({'|'.join([re.escape(s) for s in SAMPLE_FASTQS.keys()])})",  # Keep library for backward compatibility
    sample_run=f"({'|'.join([re.escape(sr) for sr in all_sample_runs])})",
    library_group=(
        f"({'|'.join([re.escape(lib) for lib in config['input']['library_groups'].keys()])})"
        if "library_groups" in config["input"]
        else ""
    ),
    resolution=f"({'|'.join([str(res) for res in config['bin']['resolutions']])})",

# --- Helper Functions ---

def argstring_to_dict(argstring):
    """
    Convert a command line argument string into a dictionary.
    
    Args:
        argstring (str): The command line argument string.
        
    Returns:
        dict: A dictionary with argument names as keys and their values.
    """
    args = shlex.split(argstring)
    keys = np.where([arg.startswith('-') for arg in args])[0]
    vals = np.where([not arg.startswith('-') for arg in args])[0]
    args_arrs = [arr for arr in np.split(args, keys) if arr.size > 0]
    argdict = {str(arr[0]): (' '.join(arr[1:]) if arr.size > 1 else True) for arr in args_arrs}
    return argdict


def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    if 'annot' not in globals():
        return []
    # Reset index to access sample_name column, then filter
    annot_reset = annot.reset_index()
    sample_runs = annot_reset[annot_reset['sample_name'] == sample_name]['sample_run'].tolist()
    return sample_runs if len(sample_runs) > 0 else []


def get_trimmed_runs_for_sample(wildcards, side=None):
    """Get all trimmed FASTQ files for a sample (all runs).
    
    Args:
        wildcards: Snakemake wildcards object
        side: Optional side (1 or 2) for paired-end reads. If None, returns both sides.
    
    Returns:
        If side is specified: list of trimmed files for that side
        If side is None: tuple of (r1_files, r2_files)
    """
    sample_runs = get_runs_for_sample(wildcards.sample)
    if side is not None:
        # Return files for specific side
        return [os.path.join(result_path, "middle_files", "trimmed", f"{sample_run}_{side}.fq.gz") 
                for sample_run in sample_runs]
    else:
        # Return both sides
        r1_files = []
        r2_files = []
        for sample_run in sample_runs:
            r1_files.append(os.path.join(result_path, "middle_files", "trimmed", f"{sample_run}_1.fq.gz"))
            r2_files.append(os.path.join(result_path, "middle_files", "trimmed", f"{sample_run}_2.fq.gz"))
        return r1_files, r2_files


def get_sample_fastqs(wildcards):
    """Return all R1/R2 FASTQ paths for a sample (all runs aggregated)."""
    return SAMPLE_FASTQS[wildcards.sample]["R1"], SAMPLE_FASTQS[wildcards.sample]["R2"]


def get_input_reads(wildcards):
    """Get input reads (trimmed or raw) for mapping."""
    if config["map"].get("trim_options"):
        return get_trimmed_runs_for_sample(wildcards)
    return get_sample_fastqs(wildcards)


def get_input_reads_for_side(wildcards):
    """Get input reads for a specific side (1 or 2) for Bowtie2."""
    if config["map"].get("trim_options"):
        return get_trimmed_runs_for_sample(wildcards, side=wildcards.side)
    fastqs = get_sample_fastqs(wildcards)
    side_idx = int(wildcards.side) - 1
    return fastqs[side_idx] if isinstance(fastqs[side_idx], list) else [fastqs[side_idx]]


def get_samtools_flags(wildcards):
    """Get samtools flags based on skip_ligation setting."""
    if SAMPLE_METADATA.get(wildcards.sample, {}).get('skip_ligation', False):
        return "-@ {threads} -bS -"
    return "-F 4 -@ {threads} -bS -"


def get_cutsite(wildcards):
    """Get cutsite for a sample."""
    return SAMPLE_METADATA.get(wildcards.sample, {}).get(
        'ligation_site', 
        config["map"].get("cutsite", "GATCGATC")
    )


def get_samples_passing_qc():
    """Get sample names that have at least one run passing QC (passqc=1).
    If passqc column doesn't exist, returns all samples (default: all pass QC)."""
    if 'annot' not in globals():
        return []
    if 'samples' not in globals():
        return []
    if 'passqc' not in annot.columns:
        return list(samples.keys())
    
    passing_samples = set()
    for sample_name in samples.keys():
        sample_runs = get_runs_for_sample(sample_name)
        for sample_run in sample_runs:
            passqc_value = annot.loc[sample_run, 'passqc']
            if passqc_value == 1 or str(passqc_value).strip() == '1':
                passing_samples.add(sample_name)
                break
    
    return list(passing_samples) if passing_samples else []


def get_units_fastqs(wildcards):
    """Get fastq files for a sample_run from PEP annotation table"""
    if 'annot' not in globals():
        raise ValueError("annot table not loaded. PEP configuration required.")
    
    u = annot.loc[wildcards.sample_run]
    
    # Support both R1/R2 and fastq1/fastq2 naming
    if 'R1' in u and 'R2' in u:
        fq1 = u["R1"]
        fq2 = u["R2"]
    elif 'fastq1' in u and 'fastq2' in u:
        fq1 = u["fastq1"]
        fq2 = u["fastq2"]
    else:
        raise ValueError(f"Sample {wildcards.sample_run} missing R1/R2 or fastq1/fastq2 columns")
    
    def expand_pep_path(path):
        """Expand PEP derive modifier paths like raw_data|filename"""
        if pd.isna(path) or not isinstance(path, str):
            return path
        if "|" in path:
            prefix, filename = path.split("|", 1)
            if prefix == "raw_data" and "paths" in config and "data_dir" in config["paths"]:
                return os.path.join(config["paths"]["data_dir"], filename)
            try:
                if 'pep_project' in globals() and 'sample_modifiers' in pep_project.config:
                    sources = pep_project.config.get('sample_modifiers', {}).get('derive', {}).get('sources', {})
                    if prefix in sources:
                        return os.path.join(sources[prefix], filename)
            except:
                pass
        if path and not os.path.isabs(path):
            return os.path.abspath(path)
        return path
    
    fq1 = expand_pep_path(fq1)
    fq2 = expand_pep_path(fq2)
    
    if pd.isna(fq2):
        fq2 = None
    
    return [fq1, fq2]


def get_all_fastqs_for_sample(sample_name):
    """Get all R1 and R2 fastq files for a sample (all runs combined)"""
    if 'annot' not in globals():
        return [], []
    
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    for sr in sample_runs:
        u = annot.loc[sr]
        if 'R1' in u:
            r1 = u["R1"]
        elif 'fastq1' in u:
            r1 = u["fastq1"]
        else:
            continue
        
        if not pd.isna(r1):
            r1_files.append(str(r1))
        
        if 'R2' in u:
            r2 = u["R2"]
        elif 'fastq2' in u:
            r2 = u["fastq2"]
        else:
            continue
        
        if not pd.isna(r2):
            r2_files.append(str(r2))
    
    return r1_files, r2_files


# --- Output-related calculations ---
min_resolution = min(config["bin"]["resolutions"])


def get_dedup_options(wildcards):
    """Get dedup options from config."""
    return config["dedup"].get("dedup_options", "")


def get_input_command(wildcards, input, threads):
    """Generate input command for merge_dedup rule."""
    return f"bgzip -dc -@ {threads-1} {input.pairs} | "


def get_phase_command(wildcards, threads):
    """Generate phase command if phasing is enabled."""
    if config.get("phase", {}).get("do_phase", False):
        return f'pairtools phase --tag-mode {config["phase"]["tag_mode"]} --phase-suffixes {" ".join(config["phase"]["suffixes"])} | pairtools sort --nproc {threads - 1} | '
    return ''


# --- Sanity Checks ---

# Some sanity checks - Phasing
if config.get('phase', {}).get('do_phase', False):
    if config['map']['mapper'] == 'chromap':
        raise ValueError(
            "Phasing is not possible with chromap, please use bwa-mem, bwa-mem2 or bwa-meme."
        )
    parse_options = argstring_to_dict(config['parse'].get('parsing_options', ''))
    if '--min-mapq' not in parse_options or parse_options['--min-mapq'] != '0':
        raise ValueError(
            "Please set '--min-mapq to 0' in the parsing options to use phasing."
        )
    if '--add-columns' not in parse_options:
        raise ValueError(
            "Please set the appropriate --add-columns argument in the parsing options to use phasing."
        )
    elif config['phase']['tag_mode'] == 'XA' and not set(['XA', 'NM', 'AS', 'XS', 'mapq']).issubset(set(parse_options['--add-columns'].split(','))):
        raise ValueError(
            "Please set '--add-columns XA,NM,AS,XS,mapq' in the parsing options to use phasing with XA tag mode."
        )
    elif config['phase']['tag_mode'] == 'XB' and not set(['XB', 'NM', 'AS', 'XS', 'mapq']).issubset(set(parse_options['--add-columns'].split(','))):
        raise ValueError(
            "Please set '--add-columns XB,NM,AS,XS,mapq' in the parsing options to use phasing with XB tag mode."
        )
    dedup_options = config["dedup"].get("dedup_options", "") 
    if '--extra-col-pair phase1 phase1' not in dedup_options or '--extra-col-pair phase2 phase2' not in dedup_options:
        raise ValueError(
            "Please add '--extra-col-pair phase1 phase2' in the dedup options to use phasing."
        )
