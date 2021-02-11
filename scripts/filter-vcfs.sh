#!/bin/bash
# FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)
#
#SBATCH --job-name=ev-filter-vcfs

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .gatk.vcf`
gatk VariantFiltration -R $GENOME -V $VCFFILE -O $root.gatk.filt.vcf \
   --filter-expression 'QD < 2.0' \
   --filter-name QD2 \
   --filter-expression 'FS > 60.0'  \
   --filter-name FS60
