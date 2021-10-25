#!/bin/bash
# Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file
# Required positional arguments:
# - file in bam format to process
# - genome file

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1
GENOME=$2

root=`basename $BAMFILE .bqsr.bam`;

picard -Xmx7g -Djava.io.tmpdir=. CollectAlignmentSummaryMetrics R=$GENOME I=$BAMFILE O=$root.stat.txt
