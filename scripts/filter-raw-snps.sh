#!/bin/bash
# filter the raw SNPs using bcftools as a template of known variable sites
# Required positional arguments:
# - raw vcf file to filter

# stop if a command fails (non-zero exit status)
set -e

VCFFILE=$1

root=`basename $VCFFILE .raw.vcf`

bcftools view -i '%QUAL>=20 && DP>5' -Oz $VCFFILE > $root.filt.vcf.gz
bcftools index $root.filt.vcf.gz
