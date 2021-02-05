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


echo "Step 1 - Index Genome Reference"
sbatch --wait scripts/index-reference-genome.sh
echo "Step 1 - Done"
echo ""


echo "Step 2 - Remove Nextera Adapters"
# create the list of fastq.gz filenames to process
ls *.fastq.gz > reads.list
./scripts/sbatch-array.sh scripts/remove-nextera-adapters.sh reads.list
echo "Step 2 - Done"
echo ""


echo "Step 3 - Map using BWA with the cleaned libraries"
# create the list of trimmed reads to process
ls *_trimmed.fq.gz > reads2.list
./scripts/sbatch-array.sh scripts/map-bwa-cleaned-libs.sh reads2.list
echo "Step 3 - Done"
echo ""


echo "Step 4 - Create BAM from SAM and make an index"
# create the list of sam files to process
ls *.sam > sams.list
./scripts/sbatch-array.sh scripts/create-bma-from-sam.sh sams.list
echo "Step 4 - Done"
echo ""


echo "Step 5 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-picard-dictionary.sh
echo "Step 5 - Done"
echo ""


echo "Pipeline - Done"

