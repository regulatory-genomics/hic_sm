import os

# Main fastp rule using PEP-based sample_run format
rule fastp:
    input:
        r1=lambda w: get_units_fastqs(w)[0],
        r2=lambda w: get_units_fastqs(w)[1],
    output:
        r1=temp(os.path.join(result_path, "middle_files", "trimmed", "{sample_run}_1.fq.gz")),
        r2=temp(os.path.join(result_path, "middle_files", "trimmed", "{sample_run}_2.fq.gz")),
        report_html=os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.html"),
        report_json=os.path.join(result_path, "Report", "fastp", "{sample_run}_fastp.json"),
    conda:
        "../envs/fastp.yml"
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/rules/fastp/{sample_run}.fastp.log"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g --thread {threads} -j {output.report_json} -h {output.report_html}"


# Aggregate trimmed files from all runs per sample
def get_trimmed_runs_for_sample(wildcards):
    """Get all trimmed FASTQ files for a sample (all runs)."""
    sample_runs = get_runs_for_sample(wildcards.sample)
    r1_files = []
    r2_files = []
    for sample_run in sample_runs:
        r1_files.append(os.path.join(result_path, "middle_files", "trimmed", f"{sample_run}_1.fq.gz"))
        r2_files.append(os.path.join(result_path, "middle_files", "trimmed", f"{sample_run}_2.fq.gz"))
    return r1_files, r2_files


rule aggregate_trimmed_fastqs:
    input:
        r1=lambda w: get_trimmed_runs_for_sample(w)[0],
        r2=lambda w: get_trimmed_runs_for_sample(w)[1],
    output:
        r1=os.path.join(result_path, "middle_files", "trimmed", "{sample}_1.fq.gz"),
        r2=os.path.join(result_path, "middle_files", "trimmed", "{sample}_2.fq.gz"),
    conda:
        "../envs/fastp.yml"
    log:
        "logs/rules/aggregate_trimmed/{sample}.log"
    threads: 4
    shell:
        r"""
        cat {input.r1} > {output.r1} 2>{log}
        cat {input.r2} > {output.r2} 2>>{log}
        """
