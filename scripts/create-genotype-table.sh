#!/bin/bash
# GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier
# to work with than a VCF)
#
#SBATCH --job-name=ev-genotype-table
#SBATCH --output=logs/ev-genotype-table-%j-%a.out

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .filt.vcf`
gatk VariantsToTable -V $VCFFILE -F POS -F TYPE -F REF -F ALT -GF AD -O $root.tab
