# Cytochrome b (cytb) Ingest Pipeline

This is the ingest pipeline for mammalian cytochrome b (cytb) sequences.

## Quick Start

```bash
# Run from the ingest directory
cd ingest
nextstrain build .
```

This will fetch human MT-CYB sequences by default. See [Configuration](#configuration) below to customize for other species or broader taxonomic coverage.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Usage

Run the ingest workflow from the ingest directory:

```bash
cd ingest
nextstrain build .
```

The workflow will create two final outputs:
- `results/metadata.tsv` - Metadata for all sequences
- `results/sequences.fasta` - The sequence data

## Configuration

The default configuration is in [`defaults/config.yaml`](./defaults/config.yaml).
The workflow is contained in [Snakefile](Snakefile) with included [rules](rules).
Each rule specifies its file inputs and output and pulls its parameters from the config.

### Important Configuration Notes

**Species-level taxonomy requirement**: The NCBI datasets tool requires species-level taxonomy IDs. The default configuration includes a test species (Homo sapiens) to demonstrate functionality. For comprehensive mammalian coverage, you have several options:

1. **Multiple species approach**: Edit `defaults/config.yaml` to include multiple species:
   ```yaml
   test_species:
     - taxon_id: "9606"   # Homo sapiens (Human)
       gene_symbol: "MT-CYB"
     - taxon_id: "10090"  # Mus musculus (Mouse)
       gene_symbol: "mt-Cytb"
     - taxon_id: "9913"   # Bos taurus (Cow)
       gene_symbol: "CYTB"
   ```

2. **Use the Entrez script**: For broader taxonomic coverage, use the provided Python script:
   ```bash
   python scripts/fetch_cytb_mammals.py \
     --email your.email@example.com \
     --max-sequences 1000 \
     --output-fasta data/sequences.fasta \
     --output-metadata data/metadata.tsv
   ```

### Customizing the pipeline

To fetch cytb sequences from specific species:
1. Edit `defaults/config.yaml`
2. Update `ncbi_taxon_id` with your species of interest
3. Update `gene_symbol` with the appropriate gene name for that species (e.g., MT-CYB, mt-Cytb, CYTB)

## Input data

### GenBank data

Cytochrome b sequences and metadata are fetched via:
- **NCBI datasets**: For species-specific queries (default method)
- **NCBI Entrez**: For broader taxonomic queries (using the provided Python script)

The pipeline searches for various cytochrome b gene name variants including:
- cytb
- cytochrome b
- MT-CYB
- CYTB
- cyt b

### Output format

The pipeline produces:
- Aligned cytochrome b sequences in FASTA format
- Metadata including organism names, taxonomy, collection dates, and geographic information

## Troubleshooting

### Common Issues

1. **"taxon requires species-level" error**: The NCBI datasets tool requires species-level taxonomy IDs. You cannot use higher-level taxa like "Mammalia" (40674). Use specific species IDs or the Entrez script for broader queries.

2. **Missing gene symbols**: Different species may use different gene symbols for cytochrome b (e.g., MT-CYB for humans, mt-Cytb for mice). Check NCBI Gene database for the correct symbol for your species.

3. **Empty results**: Ensure the gene symbol matches exactly what's in NCBI for that species. Some species may not have annotated cytochrome b genes in the datasets database.

### Alternative Approaches

If you need sequences from many mammalian species:

1. **Use the Entrez script** (`scripts/fetch_cytb_mammals.py`) which can query across all mammals
2. **Create a species list** and modify the pipeline to iterate through multiple species
3. **Download from specialized databases** like GenBank or BOLD Systems that may have broader coverage
