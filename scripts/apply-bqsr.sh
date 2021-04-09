#!/bin/bash
# APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)
#
#SBATCH --job-name=ev-apply-bqsr
#SBATCH --mem 4000
#
# Required First Argument: file containing a list of *bam2 files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bam2`
gatk --java-options -Xmx8G  ApplyBQSR -I $BAMFILE -R $GENOME --bqsr-recal-file $EVDIR/$root.table  -O $EVDIR/$root.bqsr.bam
