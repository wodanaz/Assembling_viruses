#!/bin/bash
# 'Deduplicate' or mark PCR duplicates
#
#SBATCH --job-name=ev-markdup
#SBATCH --mem 10G
#
# Required First Argument: file containing a list of *.bam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bam`

picard -Xmx7g MarkDuplicates I=$BAMFILE O=$EVDIR/$root.dedup.bam M=$EVDIR/$root.metric.txt
