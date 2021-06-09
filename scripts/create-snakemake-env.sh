#!/usr/bin/env bash
# Creates a conda environment named "snakemake" that will contain the snakemake tool

# stop if a command fails (non-zero exit status)
set -e

module load Anaconda3/2019.10-gcb02
conda create -n snakemake -c conda-forge -c bioconda mamba snakemake=6.4.1 networkx=2.5
