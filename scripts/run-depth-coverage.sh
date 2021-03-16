#!/bin/bash
# How much of the reference genome is covered by more than 1 read?
#
#SBATCH --job-name=ev-depth-cov
#SBATCH --mem 10
#
# Required First Argument: file containing a list of *bqsr.bam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bqsr.bam`;

samtools depth $BAMFILE -a > $EVDIR/$root.depth.bed
