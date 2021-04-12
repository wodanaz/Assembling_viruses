#!/bin/bash
# Compile all tab tables into one for genotype
#
#SBATCH --job-name=ev-genotype-compiler
#
# Required Environment Variables:
#  EVDIR - working directory containing *.spike.tab files

# stop if a command fails (non-zero exit status)
set -e

perl genotype_compiler.pl $EVDIR/*.spike.tab > $EVDIR/spike_genotypes.tab
cp $EVDIR/spike_genotypes.tab $EVDIR/spike_genotypes.backup.tab
sed -r 's/ /\t/g' $EVDIR/spike_genotypes.tab | sed -r 's/.spike.tab//g' | sort -k1,1n > $EVDIR/spike_genotypes.final.tab
