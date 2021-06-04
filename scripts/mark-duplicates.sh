#!/bin/bash
# 'Deduplicate' or mark PCR duplicates
# Required positional arguments:
# - file in bam to process

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1

root=`basename $BAMFILE .bam`

picard -Xmx14g -Djava.io.tmpdir=. MarkDuplicates I=$BAMFILE O=$root.dedup.bam M=$root.metric.txt
