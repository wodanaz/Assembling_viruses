#!/bin/bash
# How much of the reference genome is covered by more than 1 read?
#
#SBATCH --job-name=ev-depth-cov
#SBATCH --output=logs/ev-depth-cov-%j-%a.out
#SBATCH --mem 10

# stop if a command fails (non-zero exit status)
set -e

module load samtools/1.10-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bqsr.bam`;

samtools depth $BAMFILE -a > $root.depth.bed
