# Mammalian cytochrome b phylogeny

This is a [Nextstrain](https://nextstrain.org) build for mammalian cytochrome b sequences.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html)
for Nextstrain's suite of software tools.

## Usage

If you're unfamiliar with Nextstrain builds, you may want to follow our
[Running a Pathogen Workflow guide][] first and then come back here.

### Running the build

First, copy the ingest results to the phylogenetic data directory:

    cp ../ingest/results/sequences.fasta data/
    cp ../ingest/results/metadata.tsv data/

Then run the phylogenetic workflow:

    cd phylogenetic
    nextstrain build .

The `phylogenetic` directory will contain the workflow's intermediate files
and the final output `auspice/cytb.json`.

Once you've run the build, you can view the results with:

    nextstrain view .

## Configuration

The default configuration is in [`defaults/config.yaml`](./defaults/config.yaml).
The workflow is contained in [Snakefile](Snakefile) with included [rules](rules).
Each rule specifies its file inputs and output and pulls its parameters from the config.

### Input data

This build requires preprocessed sequence and metadata files in the `data/` directory:

* `data/sequences.fasta` - FASTA file with cytochrome b sequences
* `data/metadata.tsv` - TSV file with metadata for each sequence

These files are generated using the [ingest/](../ingest/) workflow which fetches cytochrome b sequences from GenBank.

### Configuration options

Key configuration options in `defaults/config.yaml`:

* `filter.group_by`: Groups sequences by taxonomic level (default: family)
* `filter.sequences_per_group`: Number of sequences to sample per group (default: 10)
* `filter.min_length`: Minimum sequence length (default: 1000)
* `traits.columns`: Traits to reconstruct ancestral states for (default: family)

[Nextstrain]: https://nextstrain.org
[augur]: https://docs.nextstrain.org/projects/augur/en/stable/
[auspice]: https://docs.nextstrain.org/projects/auspice/en/stable/index.html
[Installing Nextstrain guide]: https://docs.nextstrain.org/en/latest/install.html
[Running a Pathogen Workflow guide]: https://docs.nextstrain.org/en/latest/tutorials/running-a-workflow.html
