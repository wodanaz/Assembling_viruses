#!/bin/bash
# Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)
# Generates recalibration table based on various user-specified covariates
# (such as read group, reported quality score, machine cycle, and nucleotide context).
# Required positional arguments:
# - file in bam format to process
# - genome file

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1
GENOME=$2

root=`basename $BAMFILE .bam2`

tabix -p vcf $root.filt.vcf.gz -f
gatk --java-options "-Djava.io.tmpdir=." BaseRecalibrator -I $BAMFILE -R $GENOME --known-sites $root.filt.vcf.gz -O $root.table
