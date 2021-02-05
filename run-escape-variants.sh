#!/usr/bin/env bash
# Escape Variants
#
# Runs escape variants pipeline as an sbatch job
# Processes all *.fastq.gz files in the current directory
# Expects MT246667.fasta to be present

# create logs directory if necessary
mkdir -p logs

# run pipeline
sbatch scripts/escape-variants-pipeline.sh

