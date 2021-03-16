#!/bin/bash
# Remove Nextera Adapters - Array Job
#
#SBATCH --job-name=ev-rna-ary
#
# Required First Argument: file containing a list of *.fastq.gz files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE on SLURM_ARRAY_TASK_ID
READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

trim_galore --fastqc --nextera $READFILE --output_dir $EVDIR
