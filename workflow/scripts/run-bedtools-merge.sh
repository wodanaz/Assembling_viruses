#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file - Array Job
# Required positional arguments:
# - bed file to process

# stop if a command fails (non-zero exit status)
set -e

BEDFILE=$1

root=`basename $BEDFILE .depth.bed`

awk '{ if ( $3 < 5 )  print $1 "\t" $2 "\t" $2 + 1 }' ${BEDFILE} | bedtools merge |  awk '{ print $1 "\t" $2 "\t" $3 - 1  }' > $root.merged.bed
