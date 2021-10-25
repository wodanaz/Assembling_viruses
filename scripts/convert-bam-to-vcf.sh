#!/bin/bash
# BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE
# Required positional arguments:
# - bam file to convert to vcf
# - indexed genome file used by bcftools

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1
GENOME=$2

root=`basename $BAMFILE .dedup.bam`;

bcftools mpileup -Ou -f $GENOME $BAMFILE --annotate FORMAT/DPR > $root.bcf
bcftools call -vm --ploidy 1 $root.bcf > $root.raw.vcf
