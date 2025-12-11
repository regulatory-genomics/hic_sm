

def get_all_chunk_stats_for_library(wc, stat_type):
    """
    Generic input function to get all chunk-level stats for a given library.
    Handles pairstat, mapstat, and RSstat.
    """
    all_files = []
    for run in LIBRARY_RUN_FASTQS[wc.library].keys():
        chunk_ids = CHUNK_IDS[wc.library][run]

        # Define file patterns and expansion params
        if stat_type == "pairstat":
            pattern = f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}.pairstat"
            params = {"chunk_id": chunk_ids}
        elif stat_type == "mapstat":
            pattern = f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}_{{side}}.mapstat"
            params = {"chunk_id": chunk_ids, "side": ["1", "2"]}
        elif stat_type == "RSstat":
            pattern = f"{mapped_parsed_sorted_chunks_folder}/{wc.library}/{run}/{{chunk_id}}.hicpro.RSstat"
            params = {"chunk_id": chunk_ids}
        else:
            raise ValueError(f"Unknown stat_type: {stat_type}")

        all_files.extend(expand(pattern, **params))
        
    return sorted(all_files)


rule concat_stats:
    input:
        stats=lambda wc: expand(
            f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.{wc.stat_type}",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        # Use {stat_type}s to create ..._mapstats_summary.tsv, etc.
        summary=f"{pairs_library_folder}/combined_{{stat_type}}s_summary.tsv",
    log:
        "logs/concat_stats/concat_{stat_type}.log",
    wildcard_constraints:
        stat_type="mapstat|pairstat|RSstat"
    run:
        from pathlib import Path

        # This logic is now generic and works for all stat types
        library_names = [Path(f).parent.name for f in input.stats]
        input_files_str = " ".join(input.stats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )


rule concat_dedup_stats:
    input:
        expand(
            f"{pairs_library_folder}/{{library}}.{assembly}.dedup.stats",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_dedup_stats_summary.tsv",
    log:
        "logs/concat_dedup_stats.log",
    run:
        import pandas as pd
        from pathlib import Path
        
        # Dictionary to store all data
        all_stats = {}
        
        # Read each stats file
        for stats_file in input:
            library_name = Path(stats_file).stem.replace(f".{assembly}.dedup", "")
            
            # Read the two-column file
            stats_dict = {}
            with open(stats_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:  # Skip empty lines
                        parts = line.split('\t')
                        if len(parts) == 2:
                            key, value = parts
                            stats_dict[key] = value
            
            all_stats[library_name] = stats_dict
        
        # Convert to DataFrame
        df = pd.DataFrame(all_stats)
        
        # Save to TSV
        df.to_csv(output.summary, sep='\t', index=True)
        
        with open(log[0], 'w') as log_file:
            log_file.write(f"Combined {len(all_stats)} library stats files\n")
            log_file.write(f"Output written to {output.summary}\n")


rule concat_hicpro_pairstat:
    input:
        hicproPairStats=expand(
            f"{pairs_library_folder}/{{library}}.allValidPairs.mergestat",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
    output:
        summary=f"{pairs_library_folder}/combined_hicpro_pairstats_summary.tsv",
    log:
        "logs/concat_RSstat/concat_hicpro_pairstats.log",
    run:
        from pathlib import Path

        # Build library names based on the order of input files (parent folder name)
        library_names = [
            Path(f).name.replace(".allValidPairs.mergestat", "")
            for f in input.hicproPairStats
        ]
        input_files_str = " ".join(input.hicproPairStats)
        names_str = " ".join(library_names)

        shell(
            "python workflow/scripts/concat_stat.py "
            "-i {input_files_str} "
            "-n {names_str} "
            "-o {output.summary} "
            "-l {log}"
        )

rule combine_chunk_stats:
    input:
        # Pass the stat_type wildcard from the output file to the input function
        lambda wc: get_all_chunk_stats_for_library(wc, wc.stat_type)
    output:
        f"{mapped_parsed_sorted_chunks_folder}/{{library}}/assemble.{{stat_type}}"
    log:
        "logs/combine_chunk_stats/{{library}}.{stat_type}.log"
    wildcard_constraints:
        stat_type="pairstat|mapstat|RSstat"  # Ensures this rule only runs for these stats
    shell:
        "python workflow/scripts/merge_statfiles.py -f {input} > {output} 2>{log}"




rule multiqc:
    input:
        expand(
            f"{pairs_library_folder}/{{library}}.{assembly}.dedup.stats",
            library=LIBRARY_RUN_FASTQS.keys(),
        ),
        expand(
            f"{stats_library_group_folder}/{{library_group}}.{assembly}.stats",
            library_group=config["input"]["library_groups"].keys(),
        )
        if "library_groups" in config["input"] and len(config["input"]["library_groups"]) > 0
        else [],
    conda:
        "../envs/multiqc.yml"
    log:
        "logs/multiqc.log",
    params:
        input_dirs=lambda wildcards, input: list(set([Path(f).parent for f in input])),
        outdir=lambda wildcards, output: Path(output[0]).parent,
    output:
        report=f"{multiqc_folder}/multiqc_report.html",
        dir=directory(multiqc_folder),
    shell:
        r"""multiqc -f --outdir {output.dir} -m pairtools \
        {params.input_dirs} \
        >{log} 2>&1"""


