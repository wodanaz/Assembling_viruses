#!/bin/bash
# filter the raw SNPs using bcftools as a template of known variable sites
#
#SBATCH --job-name=ev-filter-snps
#SBATCH --output=logs/ev-filter-snps-%j-%a.out
#SBATCH --mem 10

# stop if a command fails (non-zero exit status)
set -e

module load bcftools/1.10.2-fasrc01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .raw.vcf`

bcftools view -i '%QUAL>=20 && DP>5' -Oz $VCFFILE > $root.filt.vcf.gz
bcftools index $root.filt.vcf.gz
