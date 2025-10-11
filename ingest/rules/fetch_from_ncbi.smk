"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = data/ncbi.ndjson

"""


# Use the Python script to fetch via Entrez
rule fetch_cytb_entrez:
    output:
        sequences=temp("data/cytb_sequences.fasta"),
        metadata=temp("data/cytb_metadata.tsv"),
    params:
        max_sequences=config.get("max_sequences", 5000),
        min_length=config.get("min_sequence_length", 1000),
        max_length=config.get("max_sequence_length", 1300),
        quality_threshold=config.get("quality_threshold", 0.95),
        email=config.get("email", "user@example.com"),
        api_key_arg=f"--api-key {config.get('api_key', '')}" if config.get("api_key", "") else "",
        refseq_only="--refseq-only" if config.get("include_refseq_only", False) else "",
        deduplicate="--deduplicate-species" if config.get("deduplicate_by_species", True) else "",
    benchmark:
        "benchmarks/fetch_cytb_entrez.txt"
    log:
        "logs/fetch_cytb_entrez.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 {workflow.basedir}/scripts/fetch_cytb_mammals.py \
            --email {params.email:q} \
            {params.api_key_arg} \
            --max-sequences {params.max_sequences:q} \
            --min-length {params.min_length:q} \
            --max-length {params.max_length:q} \
            --quality-threshold {params.quality_threshold:q} \
            {params.refseq_only} \
            {params.deduplicate} \
            --output-fasta {output.sequences:q} \
            --output-metadata {output.metadata:q}
        """


rule format_entrez_to_ndjson:
    input:
        sequences="data/cytb_sequences.fasta",
        metadata="data/cytb_metadata.tsv",
    output:
        ndjson="data/ncbi.ndjson",
    benchmark:
        "benchmarks/format_entrez_to_ndjson.txt"
    log:
        "logs/format_entrez_to_ndjson.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate passthru \
            --metadata {input.metadata:q} \
            --fasta {input.sequences:q} \
            --seq-id-column accession \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            > {output.ndjson:q}
        """
