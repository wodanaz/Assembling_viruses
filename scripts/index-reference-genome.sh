#!/usr/bin/env bash
# Index Genome Reference
# Expects GENOME environment variable to contain a path to the genome to index.
#
#SBATCH --job-name=ev-igv

# stop if a command fails (non-zero exit status)
set -e

module load bwa/0.7.12-gcb01

bwa index -a is $GENOME
