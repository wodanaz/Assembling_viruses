#!/bin/bash
# GENOTYPE TABLE (Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier
# to work with than a VCF)
#
#SBATCH --job-name=ev-genotype-table
#
# Required First Argument: file containing a list of *filt.vcf files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .filt.vcf`
gatk --java-options "-Djava.io.tmpdir=$EVDIR" VariantsToTable -V $VCFFILE -F POS -F TYPE -F REF -F ALT -GF AD -O $EVDIR/$root.tab
