#!/bin/bash
# FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)
#
#SBATCH --job-name=ev-filter-vcfs
#
# Required First Argument: file containing a list of *gatk.vcf files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .gatk.vcf`
gatk VariantFiltration -R $GENOME -V $VCFFILE -O $EVDIR/$root.gatk.filt.vcf \
   --filter-expression 'QD < 2.0' \
   --filter-name QD2 \
   --filter-expression 'FS > 60.0'  \
   --filter-name FS60
