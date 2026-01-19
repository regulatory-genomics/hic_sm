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
    HIC_TAILOR_BIN = os.path.join(workflow.basedir, "scripts", "hic-tailor")
    
    rule align_hic_tailor:
        input:
            r1=lambda wildcards: get_trimmed_runs_for_sample(wildcards)[0],
            r2=lambda wildcards: get_trimmed_runs_for_sample(wildcards)[1],
            idx=idx,
        output:
            pairs=f"{outdir}/Important_processed/Pairs/{{sample}}_hic_tailor.pairs.gz",
            bam=f"{outdir}/Important_processed/Bam/{{sample}}_hic_tailor.bam",
            flagstat=f"{outdir}/Report/HicTailor/Align/{{sample}}_hic_tailor.flagstat",
            align_log=f"{outdir}/Report/HicTailor/Align/{{sample}}_hic_tailor.json",
        params:
            idx_prefix=lambda wildcards, input: str(input.idx[0]).rsplit(".", 1)[0],
            enzyme=get_cutsite,
            r1_files=lambda wildcards, input: " ".join(input.r1) if isinstance(input.r1, list) else str(input.r1),
            r2_files=lambda wildcards, input: " ".join(input.r2) if isinstance(input.r2, list) else str(input.r2),
            fragment_bed=config["parse"]["fragment_file"],
            hic_tailor_bin=HIC_TAILOR_BIN,
        threads: 8
        resources:
            mem_mb=30000,
            runtime=20
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

            # Run hic-tailor with BAM output
            {params.hic_tailor_bin} \
                --r1 {input.r1} \
                --r2 {input.r2} \
                --index {params.idx_prefix} \
                --enzyme {params.enzyme} \
                --output {output.pairs} \
                --bam-output {output.bam} \
                --threads {threads} \
                --min-frag-len 20 \
                --fragment-bed {params.fragment_bed} \
                >{log} 2>&1

            # Generate flagstat
            samtools flagstat {output.bam} > {output.flagstat}
            
            # Copy hic-tailor log file if it exists
            if [ -f {outdir}/Important_processed/Pairs/{wildcards.sample}_hic_tailor.json ]; then
                cp {outdir}/Important_processed/Pairs/{wildcards.sample}_hic_tailor.json {output.align_log}
            else
                # Create empty log file if hic-tailor didn't create one
                touch {output.align_log}
            fi
            """

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