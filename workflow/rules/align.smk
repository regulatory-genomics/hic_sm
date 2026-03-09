# Mapping rules for BWA-mem2 and Bowtie2

# --- Index Rules ---

# Optionally build BWA-MEM2 index locally when requested.
# By default (build_bwa_index: false), we assume a shared, prebuilt index
# and treat idx files as static inputs.
if MAPPER == "bwa-mem2" and config["map"].get("build_bwa_index", False):

    rule bwaindex:
        input:
            genome=genome_path,
        output:
            idx=idx,
        params:
            bwa=MAPPER,
        conda:
            "../envs/bwa-mem2.yml"
        threads: 4
        log:
            f"logs/bwa-memx_index/{assembly}.log"
        shell:
            r"""
            # Build BWA-MEM2 index for the reference genome
            {params.bwa} index -t {threads} {input.genome} >{log} 2>&1
            """

elif MAPPER == "bowtie2":

    rule bowtie2_index:
        input:
            reference=bowtie_index_path,
        output:
            idx=idx,
        log:
            f"logs/bowtie2_index/{assembly}.log",
        params:
            extra="",
        threads: 8
        cache: True
        wrapper:
            "v3.9.0/bio/bowtie2/build"

# --- Mapping Rules ---

if MAPPER == "hic-tailor":
    import os
    HIC_TAILOR_BIN = config["map"].get(
        "hic_tailor_bin",
        os.path.join(workflow.basedir, "scripts", "hic-tailor")
    )

    rule align_hic_tailor:
        input:
            r1=lambda wildcards: get_sample_fastqs(wildcards)[0],
            r2=lambda wildcards: get_sample_fastqs(wildcards)[1],
            idx=idx,
        output:
            pairs=f"{outdir}/Important_processed/Pairs/{{sample}}.nodups.pairs.gz",
            parquet=f"{outdir}/Important_processed/Pairs/{{sample}}.hic_tailor.parquet",
            bam=f"{outdir}/Important_processed/Bam/{{sample}}_hic_tailor.bam",
            flagstat=f"{outdir}/Report/HicTailor/Align/{{sample}}_hic_tailor.flagstat",
            align_log=f"{outdir}/Report/HicTailor/Align/{{sample}}_hic_tailor.json",
        params:
            idx_prefix=lambda wildcards, input: str(input.idx[0]).rsplit(".", 1)[0],
            enzyme=get_cutsite,
            fragments=config["parse"]["fragment_file"],
            hic_tailor_bin=HIC_TAILOR_BIN,
            enable_exact_dedup=lambda wildcards: config["map"].get("hic_tailor_enable_exact_dedup", True),
            dedup_radius=lambda wildcards: config["map"].get("hic_tailor_dedup_radius", 3),
            dedup_method=lambda wildcards: config["map"].get("hic_tailor_dedup_method", "max"),
            dedup_cache_backend=lambda wildcards: config["map"].get("hic_tailor_dedup_cache_backend", "memory"),
            dedup_disk_path=lambda wildcards: config["map"].get("hic_tailor_dedup_disk_path", ""),
        threads: 16
        resources:
            mem_mb=60000,
            runtime=1500
        log:
            f"logs/hic_tailor/align/{{sample}}.log"
        benchmark:
            f"benchmarks/hic_tailor/align/{{sample}}.tsv"
        conda:
            "../envs/bwa-mem2.yml"
        shell:
            r"""
            mkdir -p {outdir}/Important_processed/Pairs
            mkdir -p {outdir}/Important_processed/Bam
            mkdir -p {outdir}/Report/HicTailor/Align

            EXTRA_DEDUP_DISK_ARGS=""
            if [ -n "{params.dedup_disk_path}" ]; then
                EXTRA_DEDUP_DISK_ARGS="--dedup-disk-path {params.dedup_disk_path}"
            fi

            EXTRA_EXACT_DEDUP_ARGS=""
            if [ "{params.enable_exact_dedup}" = "True" ] || [ "{params.enable_exact_dedup}" = "true" ]; then
                EXTRA_EXACT_DEDUP_ARGS="--enable-exact-dedup"
            fi

            {params.hic_tailor_bin} \
                -1 {input.r1} \
                -2 {input.r2} \
                --index {params.idx_prefix} \
                --enzyme '{params.enzyme}' \
                --enable-cut-mode false \
                --output {output.parquet} \
                --pairs-output {output.pairs} \
                --bam-output {output.bam} \
                --threads {threads} \
                --min-frag-len 20 \
                --fragments {params.fragments} \
                --enable-dedup \
                --dedup-radius {params.dedup_radius} \
                --dedup-method {params.dedup_method} \
                --dedup-cache-backend {params.dedup_cache_backend} \
                $EXTRA_EXACT_DEDUP_ARGS \
                $EXTRA_DEDUP_DISK_ARGS \
                >{log} 2>&1

            samtools flagstat {output.bam} > {output.flagstat}

            if [ -f {outdir}/Important_processed/Pairs/{wildcards.sample}.hic_tailor.json ]; then
                cp {outdir}/Important_processed/Pairs/{wildcards.sample}.hic_tailor.json {output.align_log}
            elif [ -f {outdir}/Important_processed/Pairs/{wildcards.sample}.json ]; then
                cp {outdir}/Important_processed/Pairs/{wildcards.sample}.json {output.align_log}
            else
                touch {output.align_log}
            fi
            """

# Legacy BWA mapping
if MAPPER == "bwa-mem2":

    rule map_chunks_bwa:
        input:
            reads=get_input_reads,
            reference=genome_path,
            idx=idx,
        params:
            bwa=MAPPER,
            extra="-SP5M",
            r1_cmd=lambda wildcards, input: " ".join(input.reads[0]) if isinstance(input.reads[0], list) else input.reads[0],
            r2_cmd=lambda wildcards, input: " ".join(input.reads[1]) if isinstance(input.reads[1], list) else input.reads[1],
        threads: 4
        output:
            f"{outdir}/Important_processed/Bam/{{sample}}.bam",
        log:
            "logs/bwa_memx/{sample}.log",
        benchmark:
            "benchmarks/bwa_memx/{sample}.tsv"
        shell:
            r"""
            # Use process substitution to concatenate multiple FASTQ files
            {params.bwa} mem {params.extra} -t {threads} {input.reference} \
            <(cat {params.r1_cmd}) \
            <(cat {params.r2_cmd}) \
            | samtools view -@ {threads} -bS - > {output} 2>{log}
            """


# Bowtie2 rescue mapping workflow is in bowtie2_rescue.smk
if MAPPER == "bowtie2":
    include: "bowtie2_rescue.smk"
