#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
# - Indexes MT246667.fasta
# - Processes all *.fastq.gz files in the current directory.
#
#SBATCH --job-name=ev-pipeline
#SBATCH --output=logs/ev-pipeline-%j.out

# stop if a command fails (non-zero exit status)
set -e

MAX_ARRAY_JOBS=4


# Run sbatch array job based on a file containing a list of filenames to process
# Called like so:
#  run_sbatch_array_job <sbatch_script_to_run> <filenames_file>
# The FILENAMES_FILE environment variable is exported to be use by the sbatch_script_to_run
run_sbatch_array_job () {
  SCRIPTNAME=$1
  # export FILENAMES_FILE for use by array job scripts
  export FILENAMES_FILE=$2

  NUM_READ_FILES=$(cat ${FILENAMES_FILE} | wc -l)
  echo "Processing ${NUM_READ_FILES} files"
  # create an array job for each file but only allow MAX_ARRAY_JOBS processes at a time
  sbatch --wait --array=1-${NUM_READ_FILES}%${MAX_ARRAY_JOBS} $SCRIPTNAME
}


echo "Pipeline - Starting"


echo "Step 1 - Index Genome Reference"
sbatch --wait scripts/index-reference-genome.sh
echo "Step 1 - Done"
echo ""


echo "Step 2 - Remove Nextera Adapters"
# create the list of fastq.gz filenames to process
ls *.fastq.gz > reads.list
run_sbatch_array_job scripts/remove-nextera-adapters.sh reads.list
echo "Step 2 - Done"
echo ""


echo "Step 3 - Map using BWA with the cleaned libraries"
# create the list of trimmed reads to process
ls *_trimmed.fq.gz > reads2.list
run_sbatch_array_job scripts/map-bwa-cleaned-libs.sh reads2.list
echo "Step 3 - Done"
echo ""


echo "Step 4 - Create BAM from SAM and make an index"
# create the list of sam files to process
ls *.sam > sams.list
run_sbatch_array_job scripts/create-bma-from-sam.sh sams.list
echo "Step 4 - Done"
echo ""


echo "Step 5 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-picard-dictionary.sh
echo "Step 5 - Done"
echo ""


echo "Pipeline - Done"

