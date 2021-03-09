#!/bin/bash
# make a table with coverage information
#
#SBATCH --job-name=ev-coverage-table
#
# Required Environment Variables:
#  EVDIR - working directory containing bams.list and bqsrs.list files

# stop if a command fails (non-zero exit status)
set -e


module load samtools

for i in `cat $EVDIR/bams.list`; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> $EVDIR/coverage.tab
done

sort -k1,1 -V $EVDIR/coverage.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > $EVDIR/coverage.raw.tab


for i in `cat $EVDIR/bqsrs.list`; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> $EVDIR/coverage.bqsr.tab
done

sort -k1,1 -V $EVDIR/coverage.bqsr.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > $EVDIR/coverage.gatk.tab
