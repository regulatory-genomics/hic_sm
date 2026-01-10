# Stitch-and-Cut Workflow
# This workflow stitches paired-end reads using FLASH2, splits at restriction sites,
# and processes through three parallel alignment streams

import os

# Get stitch_cut configuration
STITCH_CUT_CONFIG = config.get("map", {}).get("stitch_cut", {})
MIN_OVERLAP = STITCH_CUT_CONFIG.get("min_overlap", 10)
MAX_OVERLAP = STITCH_CUT_CONFIG.get("max_overlap", 150)
SAM2PAIRS_BIN = os.path.join(workflow.basedir, "scripts", "sam2pairs")
HIC_SPLIT_BIN = os.path.join(workflow.basedir, "scripts", "cut_stitched")

# Helper function to get trimmed FASTQ files for a sample (all runs aggregated)
def get_stitch_cut_input_reads(wildcards):
    """Get trimmed FASTQ files for stitch-and-cut workflow."""
    r1_files, r2_files = get_trimmed_runs_for_sample(wildcards)
    return r1_files, r2_files


# ==============================================================================
# Step 1: Stitching (FLASH2)
# ==============================================================================
rule flash2_stitch:
    input:
        r1=lambda wildcards: get_stitch_cut_input_reads(wildcards)[0],
        r2=lambda wildcards: get_stitch_cut_input_reads(wildcards)[1],
    output:
        extended=temp(f"{outdir}/Important_processed/StitchCut/{{sample}}.extendedFrags.fastq.gz"),
        not_combined_1=temp(f"{outdir}/Important_processed/StitchCut/{{sample}}.notCombined_1.fastq.gz"),
        not_combined_2=temp(f"{outdir}/Important_processed/StitchCut/{{sample}}.notCombined_2.fastq.gz"),
        hist=f"{outdir}/Report/StitchCut/{{sample}}.hist",
        histogram=f"{outdir}/Report/StitchCut/{{sample}}.histogram",
        # Save FLASH2 stdout for later parsing
        log_file=f"{outdir}/Report/StitchCut/{{sample}}.flash.log"
    params:
        min_overlap=MIN_OVERLAP,
        max_overlap=MAX_OVERLAP,
        r1_files=lambda wildcards, input: " ".join(input.r1) if isinstance(input.r1, list) else str(input.r1),
        r2_files=lambda wildcards, input: " ".join(input.r2) if isinstance(input.r2, list) else str(input.r2),
    threads: 8
    conda:
        "../envs/flash2.yml"
    log:
        f"logs/stitch_cut/flash/{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/flash/{{sample}}.tsv"
    shell:
        r"""
        # 1. Prepare inputs: concatenate multiple FASTQ files if needed (multiple runs per sample)
        cat {params.r1_files} > {wildcards.sample}_R1_combined.fastq.gz
        cat {params.r2_files} > {wildcards.sample}_R2_combined.fastq.gz

        mkdir -p {outdir}/Important_processed/StitchCut

        # 2. Run FLASH and capture its log for downstream statistics.
        #    The --output-directory flag writes files directly to the final location,
        #    so no additional mv step is needed.
        flash --threads={threads} \
            --min-overlap={params.min_overlap} \
            --max-overlap={params.max_overlap} \
            --output-prefix={wildcards.sample} \
            --output-directory={outdir}/Important_processed/StitchCut \
            --compress \
            {wildcards.sample}_R1_combined.fastq.gz \
            {wildcards.sample}_R2_combined.fastq.gz \
            > {output.log_file} 2>&1
        mv {outdir}/Important_processed/StitchCut/{wildcards.sample}.hist {output.hist}
        mv {outdir}/Important_processed/StitchCut/{wildcards.sample}.histogram {output.histogram}

        # 3. Mirror FLASH log to Snakemake log for debugging
        cat {output.log_file} > {log}

        # 4. Cleanup temporary combined input files
        rm -f {wildcards.sample}_R1_combined.fastq.gz {wildcards.sample}_R2_combined.fastq.gz
        """


# ==============================================================================
# Step 2: Cut Strategy (Separate Stream)
# ==============================================================================
rule split_stitched_reads:
    input:
        stitched=f"{outdir}/Important_processed/StitchCut/{{sample}}.extendedFrags.fastq.gz"
    output:
        cut=f"{outdir}/Important_processed/StitchCut/{{sample}}.cut.fastq.gz",
        uncut=f"{outdir}/Important_processed/StitchCut/{{sample}}.uncut.fastq.gz",
        # Stats file to record cut/uncut counts
        stats=f"{outdir}/Report/StitchCut/{{sample}}.cut_stats.txt"
    params:
        enzyme=get_cutsite,
        tool=HIC_SPLIT_BIN,
        min_len = 30,
    log:
        f"logs/stitch_cut/split/{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/split/{{sample}}.tsv"
    threads: 2
    shell:
        r"""
        # Use Rust cut_stitched binary to split stitched reads at enzyme site.
        {params.tool} \
            --input {input.stitched} \
            --cutsite {params.enzyme} \
            --min-len {params.min_len} \
            --out-cut {wildcards.sample}.cut.tmp.fastq \
            --out-uncut {wildcards.sample}.uncut.tmp.fastq \
            >{log} 2>&1

        # Count reads in temporary FASTQ files (4 lines per read)
        n_cut=$(($(wc -l < {wildcards.sample}.cut.tmp.fastq) / 4))
        n_uncut=$(($(wc -l < {wildcards.sample}.uncut.tmp.fastq) / 4))

        mkdir -p {outdir}/Report/StitchCut
        # Write simple CSV stats
        echo "Stitched_Cut_Reads,$n_cut" > {output.stats}
        echo "Stitched_Uncut_Reads,$n_uncut" >> {output.stats}

        # Compress outputs for downstream rules
        gzip -c {wildcards.sample}.cut.tmp.fastq > {output.cut}
        gzip -c {wildcards.sample}.uncut.tmp.fastq > {output.uncut}

        # Cleanup temporary files
        rm -f {wildcards.sample}.cut.tmp.fastq {wildcards.sample}.uncut.tmp.fastq
        """


# ==============================================================================
# Step 3: Alignment (Three Parallel Streams) + Stats
# ==============================================================================

# Stream A: Unstitched (Standard Pairs)
rule align_unstitched_stitchcut:
    input:
        r1=f"{outdir}/Important_processed/StitchCut/{{sample}}.notCombined_1.fastq.gz",
        r2=f"{outdir}/Important_processed/StitchCut/{{sample}}.notCombined_2.fastq.gz",
        idx=idx,
    output:
        bam=temp(f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_unstitched.bam"),
        flagstat=f"{outdir}/Report/StitchCut/Align/{{sample}}.unstitched.flagstat"
    params:
        # Stitch-and-cut uses bwa-mem2 explicitly
        bwa="bwa-mem2",
        extra="-SP5M",
        idx_prefix=lambda wildcards, input: str(input.idx[0]).rsplit(".", 1)[0],
    threads: 8
    log:
        f"logs/stitch_cut/align/unstitched_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/align/unstitched_{{sample}}.tsv"
    conda:
        "../envs/bwa-mem2.yml"
    resources:
        mem_mb= 20000,
        runtime=600
    shell:
        r"""
        {params.bwa} mem {params.extra} -t {threads} {params.idx_prefix} \
            {input.r1} {input.r2} \
            | samtools view -@ {threads} -bS - > {output.bam} 2>{log}

        mkdir -p {outdir}/Report/StitchCut/Align
        samtools flagstat {output.bam} > {output.flagstat}
        """


# Stream B: Stitched & Cut (Fragments treated as pairs)
rule align_stitched_cut:
    input:
        reads=f"{outdir}/Important_processed/StitchCut/{{sample}}.cut.fastq.gz",
        idx=idx,
    output:
        bam=temp(f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_cut.bam"),
        flagstat=f"{outdir}/Report/StitchCut/Align/{{sample}}.cut.flagstat"
    params:
        bwa="bwa-mem2",
        extra="-SP5M",
        # Derive BWA index prefix from the first index file path
        idx_prefix=lambda wildcards, input: str(input.idx[0]).rsplit(".", 1)[0],
    threads: 8
    resources:
        mem_mb= 20000,
        runtime= 600
    log:
        f"logs/stitch_cut/align/cut_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/align/cut_{{sample}}.tsv"
    conda:
        "../envs/bwa-mem2.yml"
    shell:
        r"""
        # Align as interleaved pairs (-p flag) since we split them artificially
        {params.bwa} mem {params.extra} -t {threads} -p {params.idx_prefix} \
            {input.reads} \
            | samtools view -@ {threads} -bS - > {output.bam} 2>{log}

        mkdir -p {outdir}/Report/StitchCut/Align
        samtools flagstat {output.bam} > {output.flagstat}
        """


# Stream C: Stitched & Uncut (Processed by sam2pairs)
rule align_stitched_uncut:
    input:
        reads=f"{outdir}/Important_processed/StitchCut/{{sample}}.uncut.fastq.gz",
        idx=idx,
    output:
        bam=temp(f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_uncut.bam"),
        flagstat=f"{outdir}/Report/StitchCut/Align/{{sample}}.uncut.flagstat"
    params:
        bwa="bwa-mem2",
        extra="-SP5M -T 10",
        idx_prefix=lambda wildcards, input: str(input.idx[0]).rsplit(".", 1)[0],
    threads: 8
    resources:
        mem_mb= 20000,
        runtime= 600
    log:
        f"logs/stitch_cut/align/uncut_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/align/uncut_{{sample}}.tsv"
    conda:
        "../envs/bwa-mem2.yml"
    shell:
        r"""
        # Align as single-end reads
        {params.bwa} mem {params.extra} -t {threads} {params.idx_prefix} \
            {input.reads} \
            | samtools view -@ {threads} -bS - > {output.bam} 2>{log}

        mkdir -p {outdir}/Report/StitchCut/Align
        samtools flagstat {output.bam} > {output.flagstat}
        """


# ==============================================================================
# Step 4: Parse to Pairs (Each Stream Separately - CRITICAL!)
# ==============================================================================
# IMPORTANT: Parse each BAM separately to maintain read pair adjacency.
# Merging BAMs before parsing breaks the R1/R2 pairing that pairtools expects.

# Stream A: Parse Unstitched
rule parse_unstitched_stitchcut:
    input:
        bam=f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_unstitched.bam",
        chromsizes=chromsizes_path,
    output:
        pairs=temp(f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_unstitched.pairs.gz")
    params:
        dropsam_flag="" if config["parse"].get("make_pairsam", False) else "--drop-sam",
        dropreadid_flag=(
            "--drop-readid" if config["parse"].get("drop_readid", False) else ""
        ),
        dropseq_flag="--drop-seq" if config["parse"].get("drop_seq", True) else "",
        parsing_options=config["parse"].get("parsing_options", ""),
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    log:
        f"logs/stitch_cut/parse/unstitched_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/parse/unstitched_{{sample}}.tsv"
    shell:
        r"""
        pairtools parse {params.dropsam_flag} {params.dropreadid_flag} {params.dropseq_flag} \
            {params.parsing_options} \
            -c {input.chromsizes} {input.bam} \
            | pairtools sort --nproc {threads} -o {output.pairs} \
            >{log} 2>&1
        """


# Stream B: Parse Cut (Interleaved Pairs)
rule parse_cut_stitchcut:
    input:
        bam=f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_cut.bam",
        chromsizes=chromsizes_path,
    output:
        pairs=temp(f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_cut.pairs.gz")
    params:
        dropsam_flag="" if config["parse"].get("make_pairsam", False) else "--drop-sam",
        dropreadid_flag=(
            "--drop-readid" if config["parse"].get("drop_readid", False) else ""
        ),
        dropseq_flag="--drop-seq" if config["parse"].get("drop_seq", True) else "",
        parsing_options=config["parse"].get("parsing_options", ""),
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    log:
        f"logs/stitch_cut/parse/cut_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/parse/cut_{{sample}}.tsv"
    shell:
        r"""
        pairtools parse {params.dropsam_flag} {params.dropreadid_flag} {params.dropseq_flag} \
            {params.parsing_options} \
            -c {input.chromsizes} {input.bam} \
            | pairtools sort --nproc {threads} -o {output.pairs} \
            >{log} 2>&1
        """


# Logic 2: Custom sam2pairs Parsing (Stitched-Uncut)
rule parse_uncut_sam2pair:
    input:
        bam=f"{outdir}/Important_processed/Bam/{{sample}}.stitchcut_uncut.bam",
        # ADDED: We use the unstitched file as a "Header Template"
        header_template=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_unstitched.pairs.gz"
    output:
        pairs=temp(f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_uncut.pairs.gz"),
        log_file=f"{outdir}/Report/StitchCut/{{sample}}.uncut.2pairs.log"
    params:
        prefix=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_uncut",
        stat=f"{outdir}/Report/StitchCut/{{sample}}.uncut.flash.stat",
        bin=SAM2PAIRS_BIN,
        min_ratio=STITCH_CUT_CONFIG.get("min_ratio", 0.5),
        min_mapq=STITCH_CUT_CONFIG.get("min_mapq", 10),
    threads: 4
    conda:
        "../envs/pairtools_cooler.yml"
    log:
        f"logs/stitch_cut/parse/uncut_{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/parse/uncut_{{sample}}.tsv"
    shell:
        r"""
        # 1. Steal the VALID header from the unstitched file
        # This guarantees @SQ lines and #columns match exactly for the merge step.
        zcat {input.header_template} | grep "^#" > {output.pairs}.header

        # 2. Run sam2pairs, strip headers, and ADD DUMMY MAPQ COLUMNS
        # We append "40\t40" to columns 9 and 10 to match the template header.
        {params.bin} <(samtools view -h {input.bam}) flash {params.prefix} {threads} {params.min_ratio} {params.min_mapq} 0 \
            2>{output.log_file} \
            | grep -v "^#" \
            | awk -v OFS='\t' '{{print $0, "40", "40"}}' \
            > {output.pairs}.body

        # 3. Combine Header + Body -> Sort -> Compress
        cat {output.pairs}.header {output.pairs}.body | \
        pairtools sort --nproc {threads} | bgzip -c > {output.pairs}

        # 4. Handle Logs and Stats
        if [ -f {params.prefix}.flash2pairs.log ]; then
            mv {params.prefix}.flash2pairs.log {output.log_file}
        fi

        if [ -f {params.prefix}.flash.stat ]; then
            mv {params.prefix}.flash.stat {params.stat}
        else

            echo "Warning: No stats generated" > {params.stat}
        fi

        # Cleanup
        rm -f {output.pairs}.header {output.pairs}.body
        """


# ==============================================================================
# Step 5: Merge Final Result (All Three Streams)
# ==============================================================================
rule merge_final_pairs_stitchcut:
    input:
        unstitched=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_unstitched.pairs.gz",
        cut=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_cut.pairs.gz",
        uncut=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut_uncut.pairs.gz"
    output:
        final=f"{outdir}/Important_processed/Pairs/{{sample}}.stitchcut.pairs.gz"
    threads: 8
    conda:
        "../envs/pairtools_cooler.yml"
    log:
        f"logs/stitch_cut/merge/{{sample}}.log"
    benchmark:
        f"benchmarks/stitch_cut/merge/{{sample}}.tsv"
    shell:
        r"""
        pairtools merge --nproc {threads} -o {output.final} \
            {input.unstitched} {input.cut} {input.uncut} \
            >{log} 2>&1
        """


# ==============================================================================
# Step 6: Collect Stitch-and-Cut Summary Stats
# ==============================================================================
rule collect_stitch_cut_stats:
    input:
        flash_log=f"{outdir}/Report/StitchCut/{{sample}}.flash.log",
        cut_stats=f"{outdir}/Report/StitchCut/{{sample}}.cut_stats.txt",
        fs_unstitched=f"{outdir}/Report/StitchCut/Align/{{sample}}.unstitched.flagstat",
        fs_cut=f"{outdir}/Report/StitchCut/Align/{{sample}}.cut.flagstat",
        fs_uncut=f"{outdir}/Report/StitchCut/Align/{{sample}}.uncut.flagstat"
    output:
        report=f"{outdir}/Report/{{sample}}.stitch_cut_summary.csv"
    log:
        f"logs/stitch_cut/stats/{{sample}}.log"
    run:
        import re
        import sys

        stats = {}

        try:
            # 1. Parse FLASH2 log
            with open(input.flash_log) as f:
                content = f.read()
                # Try multiple pattern variations for FLASH2 output
                total_pairs = re.search(r"Total pairs:\s+(\d+)", content) or \
                              re.search(r"Total read pairs:\s+(\d+)", content)
                combined_pairs = re.search(r"Combined pairs:\s+(\d+)", content) or \
                                 re.search(r"Successfully combined:\s+(\d+)", content)
                uncombined_pairs = re.search(r"Uncombined pairs:\s+(\d+)", content) or \
                                   re.search(r"Not combined:\s+(\d+)", content)

                stats["Total_Input_Pairs"] = int(total_pairs.group(1)) if total_pairs else 0
                stats["Stitched_Reads"] = int(combined_pairs.group(1)) if combined_pairs else 0
                stats["Unstitched_Pairs"] = int(uncombined_pairs.group(1)) if uncombined_pairs else 0

            # 2. Parse cut stats CSV
            with open(input.cut_stats) as f:
                for line in f:
                    if "," in line:
                        key, val = line.strip().split(",")
                        stats[key] = int(val)

            # 3. Parse flagstat files with error handling
            def get_mapped_reads(flagstat_file):
                try:
                    with open(flagstat_file) as fh:
                        content = fh.read()
                        m = re.search(r"(\d+) \+ \d+ mapped", content)
                        return int(m.group(1)) if m else 0
                except Exception as e:
                    print(f"Warning: Could not parse {flagstat_file}: {e}", file=sys.stderr)
                    return 0

            stats["Mapped_Unstitched"] = get_mapped_reads(input.fs_unstitched)
            stats["Mapped_Cut"] = get_mapped_reads(input.fs_cut)
            stats["Mapped_Uncut"] = get_mapped_reads(input.fs_uncut)

            # 4. Derived metrics
            total_pairs = stats["Total_Input_Pairs"]
            stats["Stitching_Rate"] = (
                stats["Stitched_Reads"] / total_pairs * 100 if total_pairs > 0 else 0
            )

            total_stitched_processed = stats.get("Stitched_Cut_Reads", 0) + stats.get(
                "Stitched_Uncut_Reads", 0
            )
            stats["Cut_Rate"] = (
                stats.get("Stitched_Cut_Reads", 0) / total_stitched_processed * 100
                if total_stitched_processed > 0
                else 0
            )

            # Alignment rates
            total_unstitched_reads = stats["Unstitched_Pairs"] * 2
            stats["Align_Rate_Unstitched"] = (
                stats["Mapped_Unstitched"] / total_unstitched_reads * 100
                if total_unstitched_reads > 0
                else 0
            )

            total_cut_frags = stats.get("Stitched_Cut_Reads", 0) * 2
            stats["Align_Rate_Cut"] = (
                stats["Mapped_Cut"] / total_cut_frags * 100
                if total_cut_frags > 0
                else 0
            )

            total_uncut_reads = stats.get("Stitched_Uncut_Reads", 0)
            stats["Align_Rate_Uncut"] = (
                stats["Mapped_Uncut"] / total_uncut_reads * 100
                if total_uncut_reads > 0
                else 0
            )

            # 5. Write summary CSV
            with open(output.report, "w") as out:
                out.write("Metric,Value\n")
                for k, v in stats.items():
                    if isinstance(v, float):
                        out.write(f"{k},{v:.2f}\n")
                    else:
                        out.write(f"{k},{v}\n")

        except Exception as e:
            print(f"Error generating stats for {wildcards.sample}: {e}", file=sys.stderr)
            # Write minimal output to prevent rule failure
            with open(output.report, "w") as out:
                out.write("Metric,Value\n")
                out.write(f"Error,{str(e)}\n")
            raise


