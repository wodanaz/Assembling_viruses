#!/bin/bash
# Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)
# Generates recalibration table based on various user-specified covariates
# (such as read group, reported quality score, machine cycle, and nucleotide context).
#
#SBATCH --job-name=ev-baserecal
#SBATCH --mem 16G
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

tabix -p vcf $EVDIR/$root.filt.vcf.gz -f
gatk --java-options "-Djava.io.tmpdir=$EVDIR" BaseRecalibrator -I $BAMFILE -R $GENOME --known-sites $EVDIR/$root.filt.vcf.gz -O $EVDIR/$root.table
