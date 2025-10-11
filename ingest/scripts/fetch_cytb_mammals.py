#!/usr/bin/env python3
"""
Fetch cytochrome b sequences from mammals using NCBI Entrez utilities.
This provides an alternative to the datasets command which requires species-level taxonomy.
"""

import argparse
import json
import sys
from Bio import Entrez, SeqIO
import time

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch cytochrome b sequences from mammals via NCBI Entrez"
    )
    parser.add_argument("--email", default="user@example.com",
        help="Email for NCBI Entrez")
    parser.add_argument("--max-sequences", type=int, default=100,
        help="Maximum number of sequences to fetch (default: 100)")
    parser.add_argument("--output-fasta", required=True,
        help="Output FASTA file path")
    parser.add_argument("--output-metadata", required=True,
        help="Output metadata TSV file path")

    return parser.parse_args()

def fetch_cytb_sequences(email, max_sequences):
    """Fetch cytochrome b sequences from mammals."""
    Entrez.email = email

    # Search query for mammalian cytochrome b
    # Includes various gene name variants
    search_query = (
        'Mammalia[Organism] AND ('
        'cytb[Gene] OR "cytochrome b"[Gene] OR "MT-CYB"[Gene] OR '
        '"CYTB"[Gene] OR "cyt b"[Gene] OR "cytochrome-b"[Gene]'
        ') AND mitochondrion[Filter]'
    )

    print(f"Searching for: {search_query}", file=sys.stderr)

    # Search for sequences
    search_handle = Entrez.esearch(
        db="nucleotide",
        term=search_query,
        retmax=max_sequences,
        idtype="acc"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_results["IdList"]
    print(f"Found {len(id_list)} sequences", file=sys.stderr)

    if not id_list:
        return []

    # Fetch sequences in batches
    sequences = []
    batch_size = 100

    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]

        print(f"Fetching batch {start//batch_size + 1}: sequences {start+1}-{end}", file=sys.stderr)

        fetch_handle = Entrez.efetch(
            db="nucleotide",
            id=batch_ids,
            rettype="gb",
            retmode="text"
        )

        for record in SeqIO.parse(fetch_handle, "genbank"):
            # Extract metadata
            metadata = {
                "accession": record.id,
                "organism": record.annotations.get("organism", ""),
                "taxonomy": ";".join(record.annotations.get("taxonomy", [])),
                "description": record.description,
                "length": len(record.seq),
                "date": record.annotations.get("date", ""),
                "country": "",
                "location": "",
            }

            # Extract location from features
            for feature in record.features:
                if feature.type == "source":
                    qualifiers = feature.qualifiers
                    metadata["country"] = qualifiers.get("country", [""])[0]
                    metadata["location"] = qualifiers.get("geo_loc_name", [""])[0]
                    break

            sequences.append((record, metadata))

        fetch_handle.close()

        # Be polite to NCBI
        time.sleep(0.34)  # Max 3 requests per second

    return sequences

def main():
    args = parse_args()

    # Fetch sequences
    sequences = fetch_cytb_sequences(args.email, args.max_sequences)

    if not sequences:
        print("No sequences found", file=sys.stderr)
        sys.exit(1)

    # Write FASTA file
    with open(args.output_fasta, "w") as fasta_file:
        for record, metadata in sequences:
            SeqIO.write(record, fasta_file, "fasta")

    # Write metadata TSV
    with open(args.output_metadata, "w") as tsv_file:
        # Write header
        headers = ["accession", "organism", "taxonomy", "description", "length", "date", "country", "location"]
        tsv_file.write("\t".join(headers) + "\n")

        # Write data
        for record, metadata in sequences:
            row = [str(metadata.get(h, "")) for h in headers]
            tsv_file.write("\t".join(row) + "\n")

    print(f"Successfully wrote {len(sequences)} sequences", file=sys.stderr)

if __name__ == "__main__":
    main()