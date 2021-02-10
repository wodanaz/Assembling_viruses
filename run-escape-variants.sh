#!/usr/bin/env bash
# Escape Variants
#
# Runs escape variants pipeline as an sbatch job
# Processes all *.fastq.gz files in the current directory
# Expects MT246667.fasta to be present

GENOME_BASE_NAME=MT246667
GENOME_INDEX_FILE="$GENOME_BASE_NAME.fasta.fai"
GENOME_DICTIONARY="$GENOME_BASE_NAME.dict"
if [[ ! -f "$GENOME_INDEX_FILE" || ! -f "$GENOME_DICTIONARY" ]]
then
    echo "ERROR: Genome index/dictionary not found."
    echo "To fix run ./setup-escape-variants.sh"
    echo ""
    exit 1
fi

# create logs directory if necessary
mkdir -p logs

# run pipeline
sbatch scripts/escape-variants-pipeline.sh

