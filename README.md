# Nextstrain repository for Cytochrome b (cytb) across Mammals

This repository contains two workflows for the analysis of mammalian cytochrome b sequences:

- [`ingest/`](./ingest) - Download cytochrome b sequences from GenBank across all mammals
- [`phylogenetic/`](./phylogenetic) - Filter sequences, align, construct phylogeny and export for visualization

Each folder contains a README.md with more information.

## Installation

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Quickstart

First, fetch cytochrome b sequences from GenBank:

    cd ingest
    nextstrain build .

Then, build the phylogenetic tree:

    cd ../phylogenetic
    cp ../ingest/results/sequences.fasta data/
    cp ../ingest/results/metadata.tsv data/
    nextstrain build .
    nextstrain view .


## Documentation

- [Running a pathogen workflow](https://docs.nextstrain.org/en/latest/tutorials/running-a-workflow.html)
- [Contributor documentation](./CONTRIBUTING.md)
