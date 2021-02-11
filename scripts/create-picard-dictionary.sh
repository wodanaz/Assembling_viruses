#!/usr/bin/env bash
# Create a dictionary file for using picard tools
# Expects GENOME environment variable to contain a path to the genome to index.
#
#SBATCH --job-name=ev-dict

# stop if a command fails (non-zero exit status)
set -e

module load picard-tools/2.4.1-gcb01
module load samtools/1.10-gcb01

GENOME_BASE_NAME=$(basename $GENOME .fasta)
GENOME_DICTIONARY="$GENOME_BASE_NAME.dict"

java -jar $PICARD_TOOLS_HOME/picard.jar  CreateSequenceDictionary R=$GENOME O=$GENOME_DICTIONARY
samtools faidx $GENOME
