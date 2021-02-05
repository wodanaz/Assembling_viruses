#!/bin/bash
# Create BAM from SAM and make an index - Array Job
#
#SBATCH --job-name=ev-mk-bam-ary
#SBATCH --output=logs/ev-mk-bam-%j-%a.out
#SBATCH --mem 1000

# stop if a command fails (non-zero exit status)
set -e

module load samtools/1.10-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $READFILE .sam`;

samtools view -Sb $READFILE | samtools sort - > $root.bam
samtools index  $root.bam
