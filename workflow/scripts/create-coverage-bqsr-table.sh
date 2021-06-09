#!/bin/bash
# make a table with coverage information
# - files in bam format to process

# stop if a command fails (non-zero exit status)
set -e

>coverage.bqsr.tab
for i in "$@"; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> coverage.bqsr.tab
done

sort -k1,1 -V coverage.bqsr.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > results/coverage.gatk.tab
