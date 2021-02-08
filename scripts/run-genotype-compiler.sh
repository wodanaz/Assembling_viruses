#!/bin/bash
# Compile all tab tables into one for genotype
#
#SBATCH --job-name=ev-genotype-compiler
#SBATCH --output=logs/ev-genotype-compiler-%j.out

# stop if a command fails (non-zero exit status)
set -e

module load perl/5.10.1-fasrc04

perl genotype_compiler.pl *.filt.tab > allgenotypes.tab
cp allgenotypes.tab allgenotypes.backup.tab
sed -r 's/ /\t/g'  allgenotypes.tab | sed -r 's/.filt.tab//g' | sort -k1,1n > allgenotypes.final.tab
