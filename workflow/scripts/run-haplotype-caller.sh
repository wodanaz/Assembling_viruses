#!/bin/bash
# HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes
# Required positional arguments:
# - bam file to process
# - genome file to use

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1
GENOME=$2

root=`basename $BAMFILE .bqsr.bam`
gatk --java-options "-Djava.io.tmpdir=." HaplotypeCaller -I $BAMFILE -R $GENOME -ploidy 1 -O $root.gatk.vcf
