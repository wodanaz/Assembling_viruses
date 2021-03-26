#!/bin/bash
# BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE
#
#SBATCH --job-name=ev-gatk-bam-to-vcf
#SBATCH --mem 10G
#
# Required First Argument: file containing a list of *dedup.bam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .dedup.bam`;

bcftools mpileup -Ou -f $GENOME $BAMFILE --annotate FORMAT/DPR > $EVDIR/$root.bcf
bcftools call -vm --ploidy 1 $EVDIR/$root.bcf > $EVDIR/$root.raw.vcf
