#!/bin/bash
# filter the raw SNPs using bcftools as a template of known variable sites
#
#SBATCH --job-name=ev-filter-snps
#SBATCH --mem 10
#
# Required First Argument: file containing a list of *raw.vcf files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .raw.vcf`

bcftools view -i '%QUAL>=20 && DP>5' -Oz $VCFFILE > $EVDIR/$root.filt.vcf.gz
bcftools index $EVDIR/$root.filt.vcf.gz
