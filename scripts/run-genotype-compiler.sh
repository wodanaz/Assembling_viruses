#!/bin/bash
# Compile all tab tables into one for genotype
#
#SBATCH --job-name=ev-genotype-compiler
#
# Required Environment Variables:
#  EVDIR - working directory containing *.filt.tab files

# stop if a command fails (non-zero exit status)
set -e

module load perl/5.10.1-fasrc04

perl genotype_compiler.pl $EVDIR/*.filt.tab > $EVDIR/allgenotypes.tab
cp $EVDIR/allgenotypes.tab $EVDIR/allgenotypes.backup.tab
sed -r 's/ /\t/g' $EVDIR/allgenotypes.tab | sed -r 's/.filt.tab//g' | sort -k1,1n > $EVDIR/allgenotypes.final.tab
