

rule parse_sort_chunks:
    input:
        bam=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.bam",
        chromsizes=chromsizes_path,
    threads: 4
    params:
        # keep_bams_command=f"| tee >(samtools view -bS > {mapped_parsed_sorted_chunks_folder}/{{sample}}.{assembly}.bam)"
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
        f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.{assembly}.pairs.gz",
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



# Since there's only one chunk per sample (no chunking), return single path
def get_pair_chunks(wildcards):
    return [f"{mapped_parsed_sorted_chunks_folder}/{wildcards.sample}.{assembly}.pairs.gz"]


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
        f"{pairs_library_folder}/{{library}}.{assembly}",
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
        f"{pairs_library_folder}/{{library}}.{assembly}",
        ".nodups.pairs.gz",
        ".unmapped.pairs.gz",
        ".dups.pairs.gz",
        ".dedup.stats",
        ".nodups.pairs.gz.px2",
    )
if config["dedup"].get("save_by_tile_dups", False):
    merge_output += [f"{pairs_library_folder}/{{library}}.{assembly}.by_tile_dups.txt"]


rule merge_dedup:
    input:
        pairs=f"{mapped_parsed_sorted_chunks_folder}/{{library}}.{assembly}.pairs.gz",
    params:
        dedup_options=lambda wildcards: config["dedup"].get("dedup_options", ""),
        max_mismatch_bp=config["dedup"]["max_mismatch_bp"],
        input_command=lambda wildcards, input, threads: (
            f"bgzip -dc -@ {threads-1} {input.pairs} | "
        ),
        phase_command=lambda wildcards, threads: (
            f'pairtools phase --tag-mode {config["phase"]["tag_mode"]} --phase-suffixes {" ".join(config["phase"]["suffixes"])} | pairtools sort --nproc {threads - 1} | '
            if config.get("phase", {}).get("do_phase", False)
            else ''
        ),
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    output:
        merge_output,
    log:
        "logs/merge_dedup/{library}.log",
    benchmark:
        "benchmarks/merge_dedup/{library}.tsv"
    shell:
        r"{params.input_command}" + r"{params.phase_command}" + dedup_command + " >{log[0]} 2>&1"


# balance_args = "--balance" if {config["bin"].get("balance", True)} else ""
# balance_options = config["bin"].get("balance_options", "")
# if balance_args and balance_options:
#     balance_args = f"{balance_args} {balance_options}"