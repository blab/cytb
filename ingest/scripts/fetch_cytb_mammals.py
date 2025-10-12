#!/usr/bin/env python3
"""
Fetch cytochrome b sequences from mammals using NCBI Entrez utilities.
This provides an alternative to the datasets command which requires species-level taxonomy.
"""

import argparse
import json
import sys
import re
import itertools
import signal
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import time
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch cytochrome b sequences from mammals via NCBI Entrez"
    )
    parser.add_argument("--email", default="user@example.com",
        help="Email for NCBI Entrez")
    parser.add_argument("--api-key", default="",
        help="Optional NCBI API key for higher rate limits")
    parser.add_argument("--max-sequences", type=int, default=100,
        help="Maximum number of sequences to fetch (default: 100)")
    parser.add_argument("--min-length", type=int, default=1000,
        help="Minimum sequence length (default: 1000)")
    parser.add_argument("--max-length", type=int, default=1300,
        help="Maximum sequence length (default: 1300)")
    parser.add_argument("--quality-threshold", type=float, default=0.95,
        help="Minimum sequence quality (max fraction of ambiguous bases, default: 0.95)")
    parser.add_argument("--refseq-only", action="store_true",
        help="Only fetch RefSeq sequences")
    parser.add_argument("--deduplicate-species", action="store_true",
        help="Keep only one sequence per species")
    parser.add_argument("--output-fasta", required=True,
        help="Output FASTA file path")
    parser.add_argument("--output-metadata", required=True,
        help="Output metadata TSV file path")

    return parser.parse_args()

def fetch_cytb_sequences(args):
    """Fetch cytochrome b sequences from mammals."""
    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # Set up signal handler for graceful interruption
    interrupted = [False]  # Use list to allow modification in nested function

    def signal_handler(signum, frame):
        print("\nInterrupted! Saving progress and exiting gracefully...", file=sys.stderr)
        interrupted[0] = True

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Enhanced search query using taxon ID for precision
    # txid40674 is the NCBI taxon ID for Mammalia
    search_query = (
        'txid40674[Orgn] AND ('
        'cytb[Gene] OR "cytochrome b"[Gene] OR "MT-CYB"[Gene] OR '
        '"CYTB"[Gene] OR "cyt b"[Gene] OR "cytochrome-b"[Gene]'
        ') AND (mitochondrion[Filter] OR "complete genome"[Title]) '
        f'AND {args.min_length}:{args.max_length}[SLEN]'
    )

    # Add RefSeq filter if requested
    if args.refseq_only:
        search_query += ' AND srcdb_refseq[PROP]'

    print(f"Searching for: {search_query}", file=sys.stderr)

    # Use History Server for efficient fetching
    search_handle = Entrez.esearch(
        db="nucleotide",
        term=search_query,
        retmax=args.max_sequences,
        usehistory="y",
        idtype="acc"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()

    count = int(search_results["Count"])
    print(f"Found {count} total sequences", file=sys.stderr)

    # Get WebEnv and QueryKey for history
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    # Limit to max_sequences
    count = min(count, args.max_sequences)
    print(f"Fetching up to {count} sequences", file=sys.stderr)

    if count == 0:
        return []

    # Fetch sequences in batches using history
    sequences = []
    batch_size = 500  # Increased batch size
    seen_species = set()  # For deduplication

    for start in range(0, count, batch_size):
        # Check for interruption
        if interrupted[0]:
            print(f"Interrupted at batch {start//batch_size + 1}. Returning {len(sequences)} sequences fetched so far.", file=sys.stderr)
            break

        end = min(start + batch_size, count)
        batch_num = start//batch_size + 1
        total_batches = (count + batch_size - 1) // batch_size
        print(f"Fetching batch {batch_num}/{total_batches}: sequences {start+1}-{end}", file=sys.stderr)

        retry_count = 0
        max_retries = 3

        while retry_count < max_retries:
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="gb",
                    retmode="text"
                )

                # Parse all records in this batch before closing handle
                batch_records = []
                try:
                    for record in SeqIO.parse(fetch_handle, "genbank"):
                        batch_records.append(record)
                except Exception as parse_error:
                    print(f"Error parsing batch {start//batch_size + 1}: {parse_error}", file=sys.stderr)
                    raise parse_error
                finally:
                    fetch_handle.close()

                break  # Success - exit retry loop

            except Exception as e:
                retry_count += 1
                print(f"Retry {retry_count}/{max_retries} for batch starting at {start}: {e}", file=sys.stderr)
                if retry_count < max_retries:
                    wait_time = 5 * retry_count  # Exponential backoff
                    print(f"Waiting {wait_time} seconds before retry...", file=sys.stderr)
                    time.sleep(wait_time)

        if retry_count == max_retries:
            print(f"Failed to fetch batch starting at {start} after {max_retries} retries", file=sys.stderr)
            print(f"Continuing with partial results...", file=sys.stderr)
            continue

        # Process the successfully fetched batch
        for record in batch_records:
            # Quality check: sequence length
            seq_len = len(record.seq)
            if seq_len < args.min_length or seq_len > args.max_length:
                continue

            # Quality check: ambiguous bases
            seq_str = str(record.seq).upper()
            n_count = seq_str.count('N') + seq_str.count('?') + seq_str.count('-')
            if n_count / seq_len > (1 - args.quality_threshold):
                continue

            organism = record.annotations.get("organism", "")

            # Deduplication by species if requested
            if args.deduplicate_species and organism in seen_species:
                continue
            seen_species.add(organism)

            # Extract enhanced metadata
            taxonomy = record.annotations.get("taxonomy", [])

            # Find order - look for known mammalian order names in taxonomy
            order = ""
            # Known mammalian orders
            mammalian_orders = {
                "Rodentia", "Primates", "Chiroptera", "Carnivora", "Artiodactyla",
                "Perissodactyla", "Lagomorpha", "Eulipotyphla", "Didelphimorphia",
                "Diprotodontia", "Monotremata", "Cingulata", "Pilosa", "Afrosoricida",
                "Macroscelidea", "Tubulidentata", "Hyracoidea", "Proboscidea",
                "Sirenia", "Scandentia", "Dermoptera", "Pholidota", "Cetacea",
                "Soricomorpha", "Erinaceomorpha", "Paucituberculata", "Microbiotheria",
                "Notoryctemorphia", "Dasyuromorphia", "Peramelemorphia", "Xenarthra"
            }

            for t in taxonomy:
                if t in mammalian_orders:
                    order = t
                    break

            metadata = {
                "accession": record.id,
                "organism": organism,
                "taxonomy": ";".join(taxonomy),
                "family_name": next((t for t in taxonomy if t.endswith("idae")), ""),
                "order": order,
                "description": record.description,
                "length": seq_len,
                "date": record.annotations.get("date", ""),
                "country": "",
                "location": "",
                "lat_lon": "",
                "specimen_voucher": "",
                "collected_by": "",
                "collection_date": "",
                "host": "",
                "isolate": "",
                "pubmed_id": "",
                "taxid": "",  # Will be populated from db_xref if available
            }

            # Extract references
            if "references" in record.annotations:
                for ref in record.annotations["references"]:
                    if hasattr(ref, 'pubmed_id'):
                        metadata["pubmed_id"] = ref.pubmed_id
                        break

            # Extract detailed metadata from features
            for feature in record.features:
                if feature.type == "source":
                    qualifiers = feature.qualifiers
                    metadata["country"] = qualifiers.get("country", [""])[0]
                    metadata["location"] = qualifiers.get("geo_loc_name", [""])[0]
                    metadata["lat_lon"] = qualifiers.get("lat_lon", [""])[0]
                    metadata["specimen_voucher"] = qualifiers.get("specimen_voucher", [""])[0]
                    metadata["collected_by"] = qualifiers.get("collected_by", [""])[0]
                    metadata["collection_date"] = qualifiers.get("collection_date", [""])[0]
                    metadata["host"] = qualifiers.get("host", [""])[0]
                    metadata["isolate"] = qualifiers.get("isolate", [""])[0]

                    # Extract taxid from db_xref (format: "taxon:12345")
                    db_xrefs = qualifiers.get("db_xref", [])
                    for xref in db_xrefs:
                        if xref.startswith("taxon:"):
                            metadata["taxid"] = xref.replace("taxon:", "")
                            break
                    break

            sequences.append((record, metadata))

        # Be polite to NCBI (with API key: 10 req/sec, without: 3 req/sec)
        if args.api_key:
            time.sleep(0.1)  # 10 requests per second with API key
        else:
            time.sleep(0.34)  # 3 requests per second without API key

        # Print progress summary every 10 batches
        if batch_num % 10 == 0:
            print(f"Progress: Fetched {len(sequences)} sequences so far ({batch_num}/{total_batches} batches)", file=sys.stderr)

    print(f"\nCompleted fetching: got {len(sequences)} total sequences", file=sys.stderr)
    return sequences


def main():
    args = parse_args()

    # Fetch sequences
    sequences = fetch_cytb_sequences(args)

    if not sequences:
        print("No sequences found", file=sys.stderr)
        sys.exit(1)

    # Write FASTA file
    with open(args.output_fasta, "w") as fasta_file:
        for record, metadata in sequences:
            SeqIO.write(record, fasta_file, "fasta")

    # Write metadata TSV
    with open(args.output_metadata, "w") as tsv_file:
        # Write header with taxid field (common_name will be added in curate pipeline)
        headers = [
            "accession", "organism", "taxid", "family_name", "order", "taxonomy",
            "description", "length", "date", "country", "location",
            "lat_lon", "specimen_voucher", "collected_by", "collection_date",
            "host", "isolate", "pubmed_id"
        ]
        tsv_file.write("\t".join(headers) + "\n")

        # Write data
        for record, metadata in sequences:

            row = [str(metadata.get(h, "")) for h in headers]
            tsv_file.write("\t".join(row) + "\n")

    # Print statistics
    print(f"\nSuccessfully wrote {len(sequences)} sequences", file=sys.stderr)

    # Species diversity statistics
    unique_species = len(set(m["organism"] for r, m in sequences))
    unique_families = len(set(m["family_name"] for r, m in sequences if m["family_name"]))
    unique_orders = len(set(m["order"] for r, m in sequences if m["order"]))

    print(f"Species diversity:", file=sys.stderr)
    print(f"  Unique species: {unique_species}", file=sys.stderr)
    print(f"  Unique families: {unique_families}", file=sys.stderr)
    print(f"  Unique orders: {unique_orders}", file=sys.stderr)

    # Geographic coverage
    countries = [m["country"] for r, m in sequences if m["country"]]
    if countries:
        unique_countries = len(set(countries))
        print(f"  Countries represented: {unique_countries}", file=sys.stderr)

if __name__ == "__main__":
    main()