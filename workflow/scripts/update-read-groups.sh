#!/bin/bash
# Add sme info for the read groups
# Required positional arguments:
# - file in bam format to process

# stop if a command fails (non-zero exit status)
set -e

BAMFILE=$1

root=`basename $BAMFILE .dedup.bam`;

picard -Xmx7g AddOrReplaceReadGroups -Djava.io.tmpdir=. I=$BAMFILE O=$root.bam2 RGSM=$root RGPU=unit1 RGLB=lib_${root} RGPL=ILLUMINA
