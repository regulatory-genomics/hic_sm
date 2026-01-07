## TODO: Combine all the statistics together
rule digest_genome:
    input:
        config["input"]["genome"]["fasta"]
    output:
        config["parse"]["fragment_file"]
    threads:1
    params:
        res_site = config["map"]["cutpattern"]
    log:
        "logs/digest_genome/digest_genome.log",
    shell:
        """
        workflow/scripts/digest_genome {input} --restriction-sites {params.res_site} \
        --out {output} > {log} 2>&1
        """

# To do : use rust version to replace
rule map2frag:
    input:
        bam=f"{outdir}/Important_processed/Bam/{{sample}}.bam",
        fragments_file = config["parse"].get("fragment_file")
    output:
        validPairs=f"{outdir}/Important_processed/Pairs/{{sample}}.hicpro.validPairs",
        validPairsRSstat = f"{outdir}/Report/Pairs/{{sample}}.hicpro.RSstat",
    threads: 4
    benchmark:
        "benchmarks/map2frag/{sample}.tsv"
    log:
        "logs/map2frag/{sample}.log",
    conda:
        "../envs/hic2frag.yml"
    resources:
        runtime= 600,
    shell:
        """
        python workflow/scripts/mapped_2hic_fragments.py \
            -f {input.fragments_file} \
            -r {input.bam} \
            -o $(dirname {output.validPairs}) \
            &> {log}
        mv "$(dirname {output.validPairs})/{wildcards.sample}.hicpro.RSstat" "{output.validPairsRSstat}" 2>/dev/null || true
        """

# Since there's only one chunk per sample, return single path
def get_all_chunks_for_sample(wildcards):
    """
    Returns the .validPairs file for the requested sample.
    """
    return [f"{outdir}/Important_processed/Pairs/{wildcards.library}.hicpro.validPairs"]
