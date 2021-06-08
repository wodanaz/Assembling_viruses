#!/bin/bash
# How much of the reference genome is covered by more than 1 read?
# Required positional arguments:
# - file in bam format to process

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1

root=`basename $BAMFILE .bqsr.bam`;

samtools depth $BAMFILE -a > $root.depth.bed
