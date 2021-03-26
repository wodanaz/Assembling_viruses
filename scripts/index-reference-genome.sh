#!/usr/bin/env bash
# Index Genome Reference
# Expects GENOME environment variable to contain a path to the genome to index.
#
#SBATCH --job-name=ev-igv

# stop if a command fails (non-zero exit status)
set -e

bwa index -a is $GENOME
