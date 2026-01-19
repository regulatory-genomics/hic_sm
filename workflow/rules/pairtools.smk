

# Only define parse_sort_chunks when NOT using hic-tailor (hic-tailor outputs pairs directly)
if MAPPER != "hic-tailor":
    rule parse_sort_chunks:
        input:
            bam=f"{outdir}/Important_processed/Bam/{{sample}}.bam",
            chromsizes=chromsizes_path,
        threads: 4
        params:
            # keep_bams_command=f"| tee >(samtools view -bS > {outdir}/Important_processed/Bam/{{sample}}.{assembly}.bam)"
            # if config["parse"]["keep_unparsed_bams"]
            # else "",
            dropsam_flag="" if config["parse"].get("make_pairsam", False) else "--drop-sam",
            dropreadid_flag=(
                "--drop-readid" if config["parse"].get("drop_readid", False) else ""
            ),
            dropseq_flag="--drop-seq" if config["parse"].get("drop_seq", True) else "",
            parsing_options=config["parse"].get("parsing_options", ""),
        conda:
            "../envs/pairtools_cooler.yml"
        output:
            f"{outdir}/Important_processed/Pairs/{{sample}}.pairs.gz",
        benchmark:
            "benchmarks/parse_sort_chunks/{sample}.tsv"
        log:
            "logs/parse_sort_chunks/{sample}.log",
        shell:
            r"""
            pairtools parse {params.dropsam_flag} {params.dropreadid_flag} {params.dropseq_flag} \
            {params.parsing_options} \
            -c {input.chromsizes} {input.bam} \
            | pairtools sort --nproc {threads} -o {output} \
            >{log[0]} 2>&1
            """

if MAPPER == "hic-tailor":
    rule parse_sort_chunks_hic_tailor:
        input:
            pairs=f"{outdir}/Important_processed/Pairs/{{sample}}_hic_tailor.pairs.gz",
            chromsizes=chromsizes_path,
        threads: 4
        conda:
            "../envs/pairtools_cooler.yml"
        output:
            f"{outdir}/Important_processed/Pairs/{{sample}}.pairs.gz",
        benchmark:
            "benchmarks/parse_sort_chunks_hic_tailor/{sample}.tsv",
        log:
            "logs/parse_sort_chunks_hic_tailor/{sample}.log",
        shell:
            r"""
            # hic-tailor may output unsorted pairs; sort + BGZIP-compress for downstream cooler/pairix
            # Use a temp file to avoid overwriting input during sorting
            pairtools sort --nproc {threads} {input.pairs} 2>>{log[0]} \
              | bgzip -@ {threads} -c > {output}.tmp
            mv {output}.tmp {output}
            """


if config["parse"]["make_pairsam"]:
    bytile_arg = (
        "--keep-parent-id --output-bytile-stats {output[7]}"
        if config["dedup"].get("save_by_tile_dups", False)
        else ""
    )
    dedup_command = (
        """pairtools dedup {params.dedup_options} \
        --max-mismatch {params.max_mismatch_bp} \
        --mark-dups \
        --output \
            >( pairtools split \
                --output-pairs {output[0]} \
                --output-sam {output[1]} \
            ) \
        --output-unmapped \
            >( pairtools split \
                --output-pairs {output[2]} \
                --output-sam {output[3]} \
            ) \
        --output-dups \
            >( pairtools split \
                --output-pairs {output[4]} \
                --output-sam {output[5]} \
            ) \
        --output-stats {output[6]} \
        """
        + bytile_arg
        + " && pairix {output[0]}"
    )
    merge_output = multiext(
        f"{outdir}/Important_processed/Pairs/{{library}}",
        ".nodups.pairs.gz",
        ".nodups.bam",
        ".unmapped.pairs.gz",
        ".unmapped.bam",
        ".dups.pairs.gz",
        ".dups.bam",
        ".dedup.stats",
        ".nodups.pairs.gz.px2",
    )
else:
    bytile_arg = (
        "--keep-parent-id --output-bytile-stats {output[5]}"
        if config["dedup"].get("save_by_tile_dups", False)
        else ""
    )
    dedup_command = (
        """pairtools dedup {params.dedup_options} \
        --max-mismatch {params.max_mismatch_bp} \
        --mark-dups \
        --output {output[0]} \
        --output-unmapped {output[1]} \
        --output-dups {output[2]} \
        --output-stats {output[3]} \
        """
        + bytile_arg
        + " >{log[0]} 2>&1 && pairix {output[0]}"
    )
    merge_output = multiext(
        f"{outdir}/Important_processed/Pairs/{{library}}",
        ".nodups.pairs.gz",
        ".unmapped.pairs.gz",
        ".dups.pairs.gz",
        ".dedup.stats",
        ".nodups.pairs.gz.px2",
    )
if config["dedup"].get("save_by_tile_dups", False):
    merge_output += [f"{outdir}/Important_processed/Pairs/{{library}}.by_tile_dups.txt"]


rule merge_dedup:
    input:
        # Primary per-library pairs; switch based on mapper and stitch_cut settings
        pairs=lambda wildcards: (
            f"{outdir}/Important_processed/Pairs/{wildcards.library}.stitchcut.pairs.gz"
            if config.get("map", {}).get("use_stitch_cut", False)
            else (
                f"{outdir}/Important_processed/Pairs/{wildcards.library}.pairs.gz"
            )
        ),
    params:
        dedup_options=get_dedup_options,
        max_mismatch_bp=config["dedup"]["max_mismatch_bp"],
        input_command=get_input_command,
        phase_command=get_phase_command,
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        merge_output,
    log:
        "logs/merge_dedup/{library}.log",
    benchmark:
        "benchmarks/merge_dedup/{library}.tsv"
    resources:
        mem_mb= 20000,
        runtime= 120,
    shell:
        r"{params.input_command}" + r"{params.phase_command}" + dedup_command + " >{log[0]} 2>&1"


rule link_dedup_stats:
    input:
        dedup_stats=f"{outdir}/Important_processed/Pairs/{{library}}.dedup.stats",
    output:
        report_stats=f"{outdir}/Report/pairtools/{{library}}.dedup.stats",
    shell:
        r"""
        mkdir -p $(dirname {output.report_stats})
        ln -sf {input.dedup_stats} {output.report_stats}
        """

