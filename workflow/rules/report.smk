
# Assemble per-side mapstat files into a single assemble.mapstat file
# Only for Bowtie2 mapper (which creates per-side mapstat files)
if config["map"]["mapper"] == "bowtie2":
    rule assemble_mapstat:
        input:
            mapstat_1=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}_1.mapstat",
            mapstat_2=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}_2.mapstat",
        output:
            assemble_mapstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.assemble.mapstat",
        shell:
            r"""
            cat {input.mapstat_1} {input.mapstat_2} > {output.assemble_mapstat}
            """

    # Copy pairstat to assemble.pairstat (for consistency)
    # Only for Bowtie2 mapper (which creates pairstat files)
    rule assemble_pairstat:
        input:
            pairstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.pairstat",
        output:
            assemble_pairstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.assemble.pairstat",
        shell:
            r"""
            cp {input.pairstat} {output.assemble_pairstat}
            """
else:
    # For non-Bowtie2 mappers, create empty assemble files
    rule assemble_mapstat:
        output:
            assemble_mapstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.assemble.mapstat",
        shell:
            r"""
            touch {output.assemble_mapstat}
            """

    rule assemble_pairstat:
        output:
            assemble_pairstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.assemble.pairstat",
        shell:
            r"""
            touch {output.assemble_pairstat}
            """

# Copy RSstat to assemble.RSstat (created by map2frag rule for all mappers)
rule assemble_RSstat:
    input:
        rsstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.hicpro.RSstat",
    output:
        assemble_rsstat=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.assemble.RSstat",
    shell:
        r"""
        cp {input.rsstat} {output.assemble_rsstat}
        """


rule concat_stats:
    input:
        stats=lambda wc: expand(
            f"{mapped_parsed_sorted_chunks_folder}/{{library}}.assemble.{wc.stat_type}",
            library=SAMPLE_FASTQS.keys(),
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
        # Extract library name from filename (e.g., "test1.assemble.mapstat" -> "test1")
        library_names = [Path(f).stem.split('.assemble')[0] for f in input.stats]
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
            library=SAMPLE_FASTQS.keys(),
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
            library=SAMPLE_FASTQS.keys(),
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


# Note: multiqc rule has been moved to qc.smk


