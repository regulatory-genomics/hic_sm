import os
from pathlib import Path

# Quality Control Rules

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


rule multiqc:
    input:
        expand(
            f"{outdir}/Important_processed/Pairs/{{library}}.dedup.stats",
            library=SAMPLE_FASTQS.keys(),
        ),
        expand(
            f"{outdir}/Report/Pairs/{{library_group}}.stats",
            library_group=config["input"]["library_groups"].keys(),
        )
        if "library_groups" in config["input"] and len(config["input"]["library_groups"]) > 0
        else [],
    conda:
        "../envs/multiqc.yml"
    log:
        "logs/multiqc.log",
    params:
        input_dirs=f"{outdir}/Report",
        outdir=f"{outdir}/Report",
    output:
        report=f"{outdir}/Report/multiqc_report.html",
    shell:
        r"""multiqc -f --outdir {params.outdir} {params.input_dirs} \
        >{log} 2>&1"""

