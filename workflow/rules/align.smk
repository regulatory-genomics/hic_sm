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

# Only define legacy BWA mapping when stitch_cut is NOT used
if MAPPER == "bwa-mem2" and not config["map"].get("use_stitch_cut", False):

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


# Only include Bowtie2 rescue workflow when stitch_cut is NOT used
if MAPPER == "bowtie2" and not config["map"].get("use_stitch_cut", False):
    # Bowtie2 rescue mapping workflow is in bowtie2_rescue.smk
    include: "bowtie2_rescue.smk"