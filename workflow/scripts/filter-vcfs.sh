#!/bin/bash
# FILTER VCFs (Filter variant calls based on INFO and/or FORMAT annotations)
# Required positional arguments:
# - vcf file to process
# - genome file to use

# stop if a command fails (non-zero exit status)
set -e

VCFFILE=$1
GENOME=$2

root=`basename $VCFFILE .gatk.vcf`
gatk --java-options "-Djava.io.tmpdir=." VariantFiltration -R $GENOME -V $VCFFILE -O results/$root.gatk.filt.vcf \
   --filter-expression 'QD < 10.0' \
   --filter-name QD2 \
   --filter-expression 'FS > 60.0'  \
   --filter-name FS60
