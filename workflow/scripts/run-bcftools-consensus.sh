#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file and trims first and last 25 bp - consensus
# Required positional arguments:
# - vcf file to process
# - genome to use

# stop if a command fails (non-zero exit status)
set -e

VCFFILE=$1
GENOME=$2

root=`basename $VCFFILE .gatk.filt.vcf`
bcftools view -Oz $VCFFILE > $VCFFILE.gz
bcftools index $VCFFILE.gz
bcftools norm -f $GENOME $VCFFILE.gz -Ob -o $root.norm.bcf
bcftools consensus -m $root.merged.bed -f $GENOME  -p ${root}_ -s ${root} -H A $VCFFILE.gz > $root.masked.fasta
seqtk trimfq -b 25 -e 25  $root.masked.fasta > $root.cleaned.fasta
