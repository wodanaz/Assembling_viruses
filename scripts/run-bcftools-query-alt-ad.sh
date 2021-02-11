#!/bin/bash
# bcftools query for "%POS %ALT" and "%POS [%AD]"
#
#SBATCH --job-name=ev-query-alt-ad

# stop if a command fails (non-zero exit status)
set -e

module load bcftools/1.10.2-fasrc01

for i in `cat vcfs2.list`; do root=`basename $i .gatk.filt.vcf`; echo $root ; bcftools query -f '%POS %ALT\n' $i | sed -r 's/ /\t/g' > $root.filt.tab ;  done
for i in `cat vcfs2.list`; do root=`basename $i .gatk.filt.vcf`; echo $root ; bcftools query -f '%POS [%AD]\n' $i | sed -r 's/,/\t/g' | awk '{ sum = $2 + $3 ; print $1 "\t" $3 / sum }' > $root.depth.tab ;done
