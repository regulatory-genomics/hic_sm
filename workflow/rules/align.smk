# Mapping rules for BWA, Bowtie2, and Chromap

# --- Index Rules ---

if MAPPER in ["bwa-mem", "bwa-mem2", "bwa-meme"]:

    rule bwaindex:
        input:
            genome=genome_path,
        output:
            idx=idx,
        params:
            bwa=MAPPER,
        threads: 1  # Only affects bwa-meme
        log:
            f"logs/bwa-memx_index/{assembly}.log",
        cache: True
        wrapper:
            "v4.6.0/bio/bwa-memx/index"

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

elif MAPPER == "chromap":

    rule chromap_index:
        input:
            genome=genome_path,
        output:
            idx=multiext(genome_path, ".chromap.index"),
        log:
            f"logs/chromap_index/{assembly}.log",
        conda:
            "../envs/chromap.yml"
        shell:
            r"chromap -i -r {input.genome} -o {output.idx} >{log} 2>&1"


# --- Mapping Rules ---

if MAPPER in ["bwa-mem", "bwa-mem2", "bwa-meme"]:

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
            f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.bam",
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


if MAPPER == "bowtie2":
    # Bowtie2 rescue mapping workflow is in bowtie2_rescue.smk
    include: "bowtie2_rescue.smk"


if MAPPER == "chromap":

    rule map_chunks_chromap:
        input:
            reads=get_input_reads,
            reference=genome_path,
            idx=multiext(genome_path, ".chromap.index"),
        params:
            extra=config["map"].get("mapping_options", ""),
            r1_files=lambda wildcards, input: ','.join(input.reads[0]) if isinstance(input.reads[0], list) else input.reads[0],
            r2_files=lambda wildcards, input: ','.join(input.reads[1]) if isinstance(input.reads[1], list) else input.reads[1],
        threads: 8
        output:
            f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.{assembly}.pairs.gz",
        log:
            "logs/chromap/{sample}.log",
        benchmark:
            "benchmarks/chromap/{sample}.tsv"
        conda:
            "../envs/chromap.yml"
        shell:
            # chromap doesn't output gzip files, so we need to pipe it to bgzip
            # It doesn't work with low memory mode, so we can't use the hic preset
            # Hence I provide all arguments manually except for --low-mem
            # Use comma-separated FASTQ files
            r"""
            chromap -e 4 -q 1 --split-alignment --pairs -x {input.idx} -r {input.reference} \
            -t {threads} {params.extra} \
            -1 {params.r1_files} \
            -2 {params.r2_files} \
            -o /dev/stdout 2>{log} | \
            bgzip > {output} \
            """
