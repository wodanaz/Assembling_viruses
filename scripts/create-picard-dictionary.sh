#!/usr/bin/env bash
# Create a dictionary file for using picard tools
# Expects GENOME environment variable to contain a path to the genome to index.
#
#SBATCH --job-name=ev-dict

# stop if a command fails (non-zero exit status)
set -e

GENOME_BASE_NAME=$(basename $GENOME .fasta)
GENOME_DICTIONARY="$GENOME_BASE_NAME.dict"

picard CreateSequenceDictionary R=$GENOME O=$GENOME_DICTIONARY
samtools faidx $GENOME
