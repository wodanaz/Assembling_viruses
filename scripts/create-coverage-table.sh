#!/bin/bash
# make a table with coverage information
#
#SBATCH --job-name=ev-coverage-table

# stop if a command fails (non-zero exit status)
set -e

module load bioawk
GENOME_LENGTH=$(bioawk -c fastx '{ print length($seq)}' $GENOME)


for i in `cat depths.list`; do
    root=`basename $i .depth.bed`
    percent=`awk -v len=$GENOME_LENGTH '{ if ( $3 == 0 )  count++ } END { print 100 - ( count*100 / len  ) }' $i `
    echo $root $percent >> table.tab
done

sort -k1,1 -V table.tab > table.sort.tab


for i in `cat bams.list`; do
    root=`basename $i .bam`
    coverage=`samtools coverage $i -H`
    echo $root $percent >> coverage.tab
done

sort -k1,1 -V coverage.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" > coverage.sort.tab


