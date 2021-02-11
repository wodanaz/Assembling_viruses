#!/bin/bash
# make a table with coverage information
#
#SBATCH --job-name=ev-coverage-table

# stop if a command fails (non-zero exit status)
set -e

for i in `cat depths.list`; do
    root=`basename $i .depth.bed`
    percent=`awk '{ if ( $3 == 0 )  count++ } END { print 100 - ( count*100 / 29871  ) }' $i `
    echo $root $percent >> table.tab
done

sort -k1,1 -V table.tab > table.sort.tab
