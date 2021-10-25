#!/usr/bin/env bash
# Create a dictionary file for using picard tools
# Required positional arguments:
# - genome file to create a picard dictionary for
# - path to create dictionary at

# stop if a command fails (non-zero exit status)
set -e

GENOME=$1
GENOME_DICTIONARY=$2

picard CreateSequenceDictionary R=$GENOME O=$GENOME_DICTIONARY
