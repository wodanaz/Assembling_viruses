#!/bin/bash
# make a table with coverage information
#
#SBATCH --job-name=ev-coverage-table

# stop if a command fails (non-zero exit status)
set -e


module load samtools

for i in `cat bams.list`; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> coverage.tab
done

sort -k1,1 -V coverage.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > coverage.raw.tab


for i in `cat bqsrs.list`; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $coverage >> coverage.bqsr.tab
done

sort -k1,1 -V coverage.bqsr.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > coverage.gatk.tab
