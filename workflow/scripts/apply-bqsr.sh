#!/bin/bash
# APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)
# Required positional arguments:
# - file in bam format to process
# - genome file

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$1
GENOME=$2

root=`basename $BAMFILE .bam2`
gatk --java-options "-Djava.io.tmpdir=." ApplyBQSR -I $BAMFILE -R $GENOME --bqsr-recal-file $root.table  -O $root.bqsr.bam
