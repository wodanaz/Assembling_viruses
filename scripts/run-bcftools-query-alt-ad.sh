#!/bin/bash
# bcftools query for "%POS %ALT" and "%POS [%AD]"
#
#SBATCH --job-name=ev-query-alt-ad
#
# Required Environment Variables:
#  EVDIR - working directory containing vcfs2.list directory

# stop if a command fails (non-zero exit status)
set -e

for i in `cat $EVDIR/vcfs2.list`; do root=`basename $i .gatk.filt.vcf`; echo $root ; bcftools query -f '%POS %ALT\n' $i | sed -r 's/ /\t/g' > $EVDIR/$root.filt.tab ;  done
for i in `cat $EVDIR/vcfs2.list`; do root=`basename $i .gatk.filt.vcf`; echo $root ; bcftools query -f '%POS [%AD]\n' $i | sed -r 's/,/\t/g' | awk '{ sum = $2 + $3 ; print $1 "\t" $3 / sum }' > $EVDIR/$root.depth.tab ;done
