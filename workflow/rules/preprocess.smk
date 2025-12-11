
rule trim:
    input:
        sample=get_raw_fastqs,
    log:
        "logs/fastp/{library}.{run}.{chunk_id}.log",
    params:
        extra=config["map"]["trim_options"],
    output:
        trimmed=[
            temp(
                f"{processed_fastqs_folder}/{{library}}/{{run}}/1.{{chunk_id}}_trimmed.fastq.gz"
            ),
            temp(
                f"{processed_fastqs_folder}/{{library}}/{{run}}/2.{{chunk_id}}_trimmed.fastq.gz"
            ),
        ],
        json=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.fastp.json",
        html=f"{mapped_parsed_sorted_chunks_folder}/{{library}}/{{run}}/{{chunk_id}}.fastp.html",
    wrapper:
        "v4.6.0/bio/fastp"
