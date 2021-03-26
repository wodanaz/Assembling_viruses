#!/bin/bash
# Create BAM from SAM and make an index - Array Job
#
#SBATCH --job-name=ev-bam-from-sam-ary
#SBATCH --mem 1000
#
# Required First Argument: file containing a list of sam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
SAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $SAMFILE .sam`;

samtools view -Sb $SAMFILE | samtools sort - > $EVDIR/$root.bam
samtools index $EVDIR/$root.bam
