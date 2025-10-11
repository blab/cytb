"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = data/ncbi.ndjson

"""


rule fetch_ncbi_dataset_package:
    params:
        ncbi_taxon_id=config["ncbi_taxon_id"],
        gene_symbol=config["gene_symbol"],
    output:
        dataset_package=temp("data/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "benchmarks/fetch_ncbi_dataset_package.txt"
    log:
        "logs/fetch_ncbi_dataset_package.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        datasets download gene symbol {params.gene_symbol:q} \
            --taxon {params.ncbi_taxon_id:q} \
            --include gene \
            --no-progressbar \
            --filename {output.dataset_package:q}
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences=temp("data/ncbi_dataset_sequences.fasta"),
    benchmark:
        "benchmarks/extract_ncbi_dataset_sequences.txt"
    log:
        "logs/extract_ncbi_dataset_sequences.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        unzip -jp {input.dataset_package:q} ncbi_dataset/data/gene.fna | \
        python3 {workflow.basedir}/scripts/fix_gene_headers.py \
        > {output.ncbi_dataset_sequences:q}
        """


rule format_ncbi_dataset_report:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv=temp("data/ncbi_dataset_report.tsv"),
    benchmark:
        "benchmarks/format_ncbi_dataset_report.txt"
    log:
        "logs/format_ncbi_dataset_report.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        # Extract accession from the annotations field in JSONL
        unzip -p {input.dataset_package:q} ncbi_dataset/data/data_report.jsonl | \
        jq -r '
            "accession\torganism\ttaxon_id\tcommon_name\tdescription\tgene_id\tsymbol\ttype\tlength" as $header |
            $header,
            (
                . |
                # Extract the genomic accession from annotations
                .annotations[0].genomicLocations[0].genomicAccessionVersion as $acc |
                "\($acc // .geneId)\t\(.taxname // "")\t\(.taxId // "")\t\(.commonName // "")\t\(.description // "")\t\(.geneId // "")\t\(.symbol // "")\t\(.type // "")\t1141"
            )
        ' > {output.ncbi_dataset_tsv:q}
        """


rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="data/ncbi_dataset_sequences.fasta",
        ncbi_dataset_tsv="data/ncbi_dataset_report.tsv",
    output:
        ndjson="data/ncbi.ndjson",
    benchmark:
        "benchmarks/format_ncbi_datasets_ndjson.txt"
    log:
        "logs/format_ncbi_datasets_ndjson.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur curate passthru \
            --metadata {input.ncbi_dataset_tsv:q} \
            --fasta {input.ncbi_dataset_sequences:q} \
            --seq-id-column accession \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            > {output.ndjson:q}
        """
