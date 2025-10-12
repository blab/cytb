"""
This part of the workflow handles curating the data into standardized
formats and expects input file

    sequences_ndjson = "data/ncbi.ndjson"

This will produce output files as

    metadata = "results/metadata.tsv"
    sequences = "results/sequences.fasta"

Parameters are expected to be defined in `config.curate`.
"""


def format_field_map(field_map: dict[str, str]) -> list[str]:
    """
    Format entries to the format expected by `augur curate --field-map`.
    When used in a Snakemake shell block, the list is automatically expanded and
    spaces are handled by quoted interpolation.
    """
    return [f'{key}={value}' for key, value in field_map.items()]


rule curate:
    input:
        sequences_ndjson="data/ncbi.ndjson",
        geolocation_rules=resolve_config_path(config["curate"]["local_geolocation_rules"]),
        annotations=resolve_config_path(config["curate"]["annotations"]),
    output:
        metadata="data/all_metadata.tsv",
        sequences="results/sequences.fasta",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        organism_regex=config["curate"]["organism_regex"],
        organism_backup_fields=config["curate"]["organism_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        genbank_location_field=config["curate"]["genbank_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        annotations_id=config["curate"]["annotations_id"],
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    benchmark:
        "benchmarks/curate.txt"
    log:
        "logs/curate.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.sequences_ndjson:q} \
            | augur curate rename \
                --field-map {params.field_map:q} \
            | augur curate normalize-strings \
            | augur curate parse-genbank-location \
                --location-field {params.genbank_location_field:q} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules:q} \
            | augur curate apply-record-annotations \
                --annotations {input.annotations:q} \
                --id-field {params.annotations_id:q} \
                --output-metadata {output.metadata:q} \
                --output-fasta {output.sequences:q} \
                --output-id-field {params.id_field:q} \
                --output-seq-field {params.sequence_field:q}
        """


rule fetch_common_names:
    """Fetch common names for species using NCBI Taxonomy database"""
    input:
        metadata = "data/all_metadata.tsv",
        script = "scripts/fetch_common_names.py"
    output:
        metadata = temp("data/all_metadata_common_names.tsv")
    params:
        email = config.get("email", ""),
        api_key = config.get("api_key", "")
    benchmark:
        "benchmarks/fetch_common_names.txt"
    log:
        "logs/fetch_common_names.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 {input.script:q} \
            --email {params.email:q} \
            $(if [ -n "{params.api_key}" ]; then echo "--api-key {params.api_key:q}"; fi) \
        < {input.metadata:q} \
        > {output.metadata:q}
        """


rule add_metadata_columns:
    """Add columns to metadata
    Notable columns:
    - [NEW] url: URL linking to the NCBI GenBank record ('https://www.ncbi.nlm.nih.gov/nuccore/*').
    """
    input:
        metadata = "data/all_metadata_common_names.tsv"
    output:
        metadata = temp("data/all_metadata_added.tsv")
    params:
        accession=config['curate']['genbank_accession']
    benchmark:
        "benchmarks/add_metadata_columns.txt"
    log:
        "logs/add_metadata_columns.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk mutate2 -t \
          -n url \
          -e '"https://www.ncbi.nlm.nih.gov/nuccore/" + ${params.accession:q}' \
          {input.metadata:q} \
        > {output.metadata:q}
        """


rule subset_metadata:
    input:
        metadata="data/all_metadata_added.tsv",
    output:
        subset_metadata="results/metadata.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    benchmark:
        "benchmarks/subset_metadata.txt"
    log:
        "logs/subset_metadata.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.metadata_fields:q} \
            {input.metadata:q} > {output.subset_metadata:q}
        """
