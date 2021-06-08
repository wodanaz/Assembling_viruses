#!/bin/bash
# make a table with coverage information
# Required positional arguments:
# - files in bam format to process

# stop if a command fails (non-zero exit status)
set -e

>coverage.tab
for i in "$@"; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> coverage.tab
done

sort -k1,1 -V coverage.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > coverage.raw.tab
