#!/bin/bash
# Compile all tab tables into one for genotype
# - tab files to process
# Required environment variable PERLSCRIPT that points to depth compiler

# stop if a command fails (non-zero exit status)
set -e

perl $PERLSCRIPT "$@" > spike_genotypes.tab
cp spike_genotypes.tab spike_genotypes.backup.tab
sed -r 's/ /\t/g' spike_genotypes.tab | sed -r 's/.spike.tab//g' | sort -k1,1n > results/spike_genotypes.final.tab
