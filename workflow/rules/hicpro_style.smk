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
        bam=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.bam",
        fragments_file = config["parse"].get("fragment_file")
    output:
        validPairs=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.hicpro.validPairs",
        validPairsRSstat = f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.hicpro.RSstat",
    threads: 4
    benchmark:
        "benchmarks/map2frag/{library}.{run}.{chunk_id}.tsv"
    log:
        "logs/map2frag/{library}.{run}.{chunk_id}.log",
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

# Find out the chunk ids for each library and run - since we don't know them beforehand
def get_all_chunks_for_library(wildcards):
    """
    Gathers all chunk-level .validPairs files from ALL runs that
    belong to the requested library.
    """
    all_chunk_paths = []
    # Get all run keys (e.g., 'run1', 'run2') for the given library
    runs_for_library = LIBRARY_RUN_FASTQS[wildcards.library].keys()

    # Loop over each run to find its chunks
    for run in runs_for_library:
        chunk_ids = CHUNK_IDS[wildcards.library][run]

        # Generate the paths for this run's chunks and add them to the main list
        run_chunk_paths = expand(
            f"{mapped_parsed_sorted_chunks_folder}/{wildcards.library}/{run}/{{chunk_id}}.hicpro.validPairs",
            library=wildcards.library,
            run=run,
            chunk_id=chunk_ids
        )
        all_chunk_paths.extend(run_chunk_paths)

    return all_chunk_paths

rule merge_library_pairs:
    input:
        get_all_chunks_for_library
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
