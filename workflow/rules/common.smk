# Copied from https://github.com/mirnylab/distiller-sm/blob/master/_distiller_common.py
import os
import numpy as np
import pandas as pd
import shlex

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



# Helper functions for PEP-based workflows
def get_runs_for_sample(sample_name):
    """Get all runs for a given sample name"""
    if 'annot' not in globals():
        return []
    # Reset index to access sample_name column, then filter
    annot_reset = annot.reset_index()
    sample_runs = annot_reset[annot_reset['sample_name'] == sample_name]['sample_run'].tolist()
    return sample_runs if len(sample_runs) > 0 else []


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
