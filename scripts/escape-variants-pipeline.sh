#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
# - Indexes MT246667.fasta
# - Processes all *.fastq.gz files in the current directory.
#
#SBATCH --job-name=ev-pipeline
#SBATCH --output=logs/ev-pipeline-%j.out

# stop if a command fails (non-zero exit status)
set -e

echo "Pipeline - Starting"
MAX_ARRAY_JOBS=4


echo "Step 1 - Index Genome Reference"

sbatch --wait scripts/index-reference-genome.sh

echo "Step 1 - Done"


echo "Step 2 - Remove Nextera Adapters"

# create the list of fastq.gz files to process
ls *.fastq.gz > reads.list

# count how many fastq.gz files we have
NUM_READ_FILES=$(cat reads.list | wc -l)

echo "Removing Nextera Adapaters from $NUM_READ_FILES files"
# create an array job for each fastq.gz file but only allow MAX_ARRAY_JOBS processes at a time
sbatch --wait --array=1-${NUM_READ_FILES}%${MAX_ARRAY_JOBS} scripts/remove-nextera-adapters.sh

echo "Step 2 - Done"


echo "Pipeline - Done"

