rule calculate_complexity:
    input:
        dedup_stats=f"{outdir}/Important_processed/Pairs/{{library}}.dedup.stats",
    output:
        complexity=f"{outdir}/Report/pairtools/{{library}}.complexity.tsv",
    params:
        optical_rate_arg=lambda wildcards: (
            f"--optical-rate {config.get('complexity', {}).get('optical_rate', '')}"
            if config.get("complexity", {}).get("optical_rate") is not None
            else ""
        ),
    log:
        "logs/calculate_complexity/{library}.log",
    benchmark:
        "benchmarks/calculate_complexity/{library}.tsv"
    shell:
        r"""
        python workflow/scripts/calculate_complexity.py \
            -i {input.dedup_stats} \
            -o {output.complexity} \
            --library {wildcards.library} \
            {params.optical_rate_arg} \
            >{log[0]} 2>&1
        """


rule dist_vs_contact:
    input:
        mcool=f"{outdir}/Important_processed/Cooler/{{library}}.{{filter_name}}.{{min_resolution}}.mcool",
    output:
        dist_stat = f"{outdir}/downstream_res/downstream_dist/{{library}}.{{filter_name}}.{{min_resolution}}_loglog_fits.csv",
        pic1 = f"{outdir}/downstream_res/downstream_dist/{{library}}.{{filter_name}}.{{min_resolution}}_hexbin_all_arms.png",
        pic2 = f"{outdir}/downstream_res/downstream_dist/{{library}}.{{filter_name}}.{{min_resolution}}_per_arm.png",
        report_dist_stat=f"{outdir}/Report/downstream/{{library}}.{{filter_name}}.{{min_resolution}}_loglog_fits.csv",
    log:
        "logs/downstream_dist/{library}.{filter_name}.{min_resolution}.log",
    benchmark:
        "benchmarks/downstream_dist/{library}.{filter_name}.{min_resolution}.tsv"
    conda:
        "../envs/downstream_plot.yml"
    params:
        genome = config["input"]["genome"]["assembly_name"],
        out_dir = f"{outdir}/downstream_res/downstream_dist"
    threads: 8
    shell:
        r"""
        python workflow/scripts/calculate_dist_contact.py -i {input.mcool} \
        -g {params.genome} -r {wildcards.min_resolution} -t {threads} \
        -o {params.out_dir} \
        -p {wildcards.library}.{wildcards.filter_name}.{wildcards.min_resolution} >{log[0]} 2>&1
        mkdir -p $(dirname {output.report_dist_stat})
        ln -sf {output.dist_stat} {output.report_dist_stat}
        """

rule mustache_loop_detection:
    input:
        mcool=f"{outdir}/Important_processed/Cooler/{{library}}.{{filter_name}}.{min_resolution}.mcool",
    output:
        loops = f"{outdir}/Important_processed/Loops/{{library}}.{{filter_name}}.{{resolution}}_loops.tsv"
    log:
        "logs/downstream_loops/{library}.{filter_name}.{resolution}.log",
    benchmark:
        "benchmarks/downstream_loops/{library}.{filter_name}.{resolution}.tsv"
    conda:
        "../envs/mustache.yml"
    params:
        threshold = lambda wildcards: config["mustache"]["thresholds"].get(wildcards.resolution, "0.1"),
        sparse = lambda wildcards: config["mustache"]["sparse_params"].get(wildcards.resolution, "0.88"),
        mustache_path = "workflow/scripts/mustache.py",
        prefix = "{library}.{filter_name}.{resolution}",
        out_dir = f"{outdir}/Important_processed/Loops"
    threads: 5
    resources:
        mem_mb= 100000,
        time= "6:00:00"
    shell:
        r"""
        bash workflow/scripts/mustache_loop_detection.sh \
        -i {input.mcool} \
        -r {wildcards.resolution} \
        -o {params.out_dir} \
        -p {params.prefix} \
        -t {threads} \
        -pt {params.threshold} \
        -st {params.sparse} \
        -m {params.mustache_path} >{log[0]} 2>&1
        """


rule summarize_loop_counts:
    input:
        lambda wildcards: expand(
            f"{outdir}/Important_processed/Loops/{{library}}.{{filter_name}}.{{resolution}}_loops.tsv",
            library=wildcards.library,
            filter_name=wildcards.filter_name,
            resolution=config["bin"]["resolution_for_mustache"],
        ),
    output:
        loop_counts=f"{outdir}/Report/downstream/{{library}}.{{filter_name}}.loop_counts.tsv",
    log:
        "logs/summarize_loop_counts/{library}.{filter_name}.log",
    benchmark:
        "benchmarks/summarize_loop_counts/{library}.{filter_name}.tsv"
    shell:
        r"""
        python workflow/scripts/count_loops.py \
            -i {input} \
            -o {output.loop_counts} \
            --library {wildcards.library} \
            --filter {wildcards.filter_name} \
            >{log[0]} 2>&1
        """