#!/usr/bin/env python3
"""
Fetch common names for species using NCBI Taxonomy database.
Reads metadata TSV from stdin, adds common_name column, outputs to stdout.

Usage as augur curate tool:
    cat metadata.tsv | python3 fetch_common_names.py --email user@example.com > metadata_with_common_names.tsv
"""

import sys
import time
import argparse
import pandas as pd
from io import StringIO
from Bio import Entrez


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add common names to metadata using NCBI Taxonomy",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api-key", help="NCBI API key for higher rate limits")
    parser.add_argument("--batch-size", type=int, default=800,
                       help="Batch size for bulk taxonomy fetches")
    parser.add_argument("--taxid-field", default="taxid",
                       help="Name of the column containing NCBI taxids")
    parser.add_argument("--organism-field", default="organism_name",
                       help="Name of the column containing scientific names")
    parser.add_argument("--output-field", default="common_name",
                       help="Name of the output column for common names")
    return parser.parse_args()


def fetch_common_names_bulk(taxids, email, api_key=None, batch_size=800):
    """
    Fetch common names for a list of taxonomy IDs using bulk efetch.
    Returns dict mapping taxid -> common_name.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    if not taxids:
        return {}

    # Remove duplicates while preserving order
    unique_taxids = list(dict.fromkeys(taxids))

    print(f"Fetching common names for {len(unique_taxids)} unique TaxIDs in bulk...",
          file=sys.stderr)

    common_by_taxid = {}

    for i in range(0, len(unique_taxids), batch_size):
        batch = unique_taxids[i:i+batch_size]
        batch_end = min(i+batch_size, len(unique_taxids))
        print(f"  Processing batch {i//batch_size + 1}: TaxIDs {i+1}-{batch_end}",
              file=sys.stderr)

        try:
            fh = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            recs = Entrez.read(fh)
            fh.close()

            if isinstance(recs, list):
                for rec in recs:
                    tid = str(rec.get("TaxId"))
                    other = rec.get("OtherNames") or {}
                    common = None

                    # Prefer GenBank's standardized common name
                    if "GenbankCommonName" in other:
                        common = other["GenbankCommonName"]
                    elif "CommonName" in other:
                        cn = other["CommonName"]
                        if isinstance(cn, list) and cn:
                            common = cn[0]
                        elif cn and not isinstance(cn, list):
                            common = cn

                    if common:
                        common_by_taxid[tid] = common

        except Exception as e:
            print(f"  Warning: taxonomy efetch batch failed ({e})", file=sys.stderr)

        # Rate-limit: 10 req/s with key; otherwise ~3 req/s
        time.sleep(0.1 if api_key else 0.34)

    print(f"  Found common names for {len(common_by_taxid)} TaxIDs", file=sys.stderr)
    return common_by_taxid


def resolve_missing_taxids(scientific_names, email, api_key=None):
    """
    For species without TaxIDs, try to resolve them via scientific name.
    Returns dict mapping scientific_name -> taxid.
    """
    if not scientific_names:
        return {}

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    unique_names = list(dict.fromkeys(scientific_names))
    print(f"Resolving TaxIDs for {len(unique_names)} species via name lookup...",
          file=sys.stderr)

    name_to_taxid = {}

    # Process in batches of 200 using OR queries
    for i in range(0, len(unique_names), 200):
        chunk = unique_names[i:i+200]
        query = " OR ".join(f'"{name}"[Scientific Name]' for name in chunk)

        try:
            handle = Entrez.esearch(db="taxonomy", term=query, retmax=100000)
            result = Entrez.read(handle)
            handle.close()

            ids = result.get("IdList", [])
            if ids:
                # Get summaries to map IDs back to names
                handle_summary = Entrez.esummary(db="taxonomy", id=",".join(ids), retmode="xml")
                summaries = Entrez.read(handle_summary)
                handle_summary.close()

                for rec in summaries:
                    sci_name = rec.get("ScientificName")
                    taxid = rec.get("TaxId")
                    if sci_name and taxid and sci_name in unique_names:
                        name_to_taxid[sci_name] = str(taxid)

            # Rate limiting
            time.sleep(0.1 if api_key else 0.34)

        except Exception as e:
            print(f"  Warning: name->taxid batch failed ({e})", file=sys.stderr)

    print(f"  Resolved {len(name_to_taxid)} TaxIDs from scientific names", file=sys.stderr)
    return name_to_taxid


def main():
    args = parse_args()

    # Read metadata from stdin
    print("Reading metadata from stdin...", file=sys.stderr)
    input_data = sys.stdin.read()
    df = pd.read_csv(StringIO(input_data), sep='\t')

    # Check for required columns
    if args.taxid_field not in df.columns:
        print(f"Warning: taxid field '{args.taxid_field}' not found in metadata",
              file=sys.stderr)
        # Just pass through without modification
        df.to_csv(sys.stdout, sep='\t', index=False)
        return

    if args.organism_field not in df.columns:
        print(f"Error: organism field '{args.organism_field}' not found in metadata",
              file=sys.stderr)
        sys.exit(1)

    # Collect TaxIDs and track which rows need lookup
    taxids_to_fetch = []
    missing_taxid_names = []

    for idx, row in df.iterrows():
        taxid = row[args.taxid_field]
        organism = row[args.organism_field]

        # Check if we have a valid taxid
        if pd.notna(taxid) and str(taxid).strip() and str(taxid).strip() != '':
            taxids_to_fetch.append(str(taxid))
        elif pd.notna(organism) and str(organism).strip():
            missing_taxid_names.append(str(organism))

    # First, resolve missing TaxIDs from scientific names
    name_to_taxid = {}
    if missing_taxid_names:
        name_to_taxid = resolve_missing_taxids(missing_taxid_names, args.email, args.api_key)
        # Add resolved TaxIDs to fetch list
        taxids_to_fetch.extend(name_to_taxid.values())

    # Fetch common names for all TaxIDs
    common_by_taxid = {}
    if taxids_to_fetch:
        common_by_taxid = fetch_common_names_bulk(taxids_to_fetch, args.email,
                                                  args.api_key, args.batch_size)

    # Add common_name column to dataframe
    common_names = []
    for idx, row in df.iterrows():
        taxid = str(row[args.taxid_field]) if pd.notna(row[args.taxid_field]) else None
        organism = str(row[args.organism_field]) if pd.notna(row[args.organism_field]) else None

        common_name = None

        # First try direct taxid lookup
        if taxid and taxid in common_by_taxid:
            common_name = common_by_taxid[taxid]
        # Then try resolved taxid from scientific name
        elif organism and organism in name_to_taxid:
            resolved_taxid = name_to_taxid[organism]
            if resolved_taxid in common_by_taxid:
                common_name = common_by_taxid[resolved_taxid]

        # Fall back to scientific name if no common name found
        if not common_name:
            common_name = organism if organism else ""

        common_names.append(common_name)

    # Add the column
    df[args.output_field] = common_names

    # Output to stdout
    print(f"Writing metadata with {args.output_field} column to stdout...", file=sys.stderr)
    df.to_csv(sys.stdout, sep='\t', index=False)

    # Print statistics
    n_with_common = sum(1 for cn, org in zip(common_names, df[args.organism_field])
                       if cn and cn != org)
    print(f"Added common names for {n_with_common}/{len(df)} records", file=sys.stderr)


if __name__ == "__main__":
    main()