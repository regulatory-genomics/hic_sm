
import os

result_path = config.get("output_base_dir", "results")


rule fastp:
    input:
        r1=lambda w: get_units_fastqs(w)[0],
        r2=lambda w: get_units_fastqs(w)[1],
    output:
        r1=temp(os.path.join(result_path, "middle_files", "trimmed", "{sample_run}_1.fq.gz")),
        r2=temp(os.path.join(result_path, "middle_files", "trimmed", "{sample_run}_2.fq.gz")),
        report_html=os.path.join(result_path, "report", "fastp", "{sample_run}_fastp.html"),
        report_json=os.path.join(result_path, "report", "fastp", "{sample_run}_fastp.json"),
    conda:
        "../envs/fastp.yaml"
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/rules/fastp/{sample_run}.fastp.json"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --detect_adapter_for_pe --trim_poly_g  --thread {threads} -j {output.report_json} -h {output.report_html}"

