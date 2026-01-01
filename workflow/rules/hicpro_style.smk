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
        bam=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.bam",
        fragments_file = config["parse"].get("fragment_file")
    output:
        validPairs=f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.hicpro.validPairs",
        validPairsRSstat = f"{mapped_parsed_sorted_chunks_folder}/{{sample}}.hicpro.RSstat",
    threads: 4
    benchmark:
        "benchmarks/map2frag/{sample}.tsv"
    log:
        "logs/map2frag/{sample}.log",
    conda:
        "../envs/hic2frag.yml"
    shell:
        """
        python workflow/scripts/mapped_2hic_fragments.py \
            -f {input.fragments_file} \
            -r {input.bam} \
            -o $(dirname {output.validPairs}) \
            &> {log}
        """

# Since there's only one chunk per sample, return single path
def get_all_chunks_for_sample(wildcards):
    """
    Returns the .validPairs file for the requested sample.
    """
    return [f"{mapped_parsed_sorted_chunks_folder}/{wildcards.library}.hicpro.validPairs"]

rule merge_library_pairs:
    input:
        get_all_chunks_for_sample
    threads: 2
    output:
        f"{pairs_library_folder}/{{library}}.allValidPairs",
        f"{pairs_library_folder}/{{library}}.allValidPairs.mergestat",
    log:
        "logs/merge_validpair/{library}.log",
    benchmark:
        "benchmarks/merge_validpair/{library}.tsv"
    params:
        prefix=f"{pairs_library_folder}/{{library}}"
    shell:
        "bash workflow/scripts/hicpro_merge_validpairs.sh -d -p {params.prefix} {input} &> {log}"
