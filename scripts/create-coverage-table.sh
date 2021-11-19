#!/bin/bash
# make a table with coverage information for a bam file
# Required positional arguments:
# Required positional arguments:
# - input file in bam format to process
# - output filename

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1
OUTFILE=$2

root=`basename $BAMFILE .bam`
coverage=`samtools coverage $BAMFILE -H`
echo $root $coverage > $OUTFILE
