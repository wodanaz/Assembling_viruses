#!/bin/bash
# Create BAM from SAM and make an index
# Required positional arguments:
# - file in sam format used to create a bam file

# stop if a command fails (non-zero exit status)
set -e

SAMFILE=$1
OUTFILE=$2

samtools view -Sb $SAMFILE | samtools sort - > $OUTFILE
samtools index $OUTFILE
