#!/bin/bash
# Compile all tab tables into one for depth
#
#SBATCH --job-name=ev-depth-compiler
#SBATCH --output=logs/ev-depth-compiler-%j.out

# stop if a command fails (non-zero exit status)
set -e

module load perl/5.10.1-fasrc04

perl depth_compiler.pl *.depth.tab > alldepths.tab
cp alldepths.tab alldepths.backup.tab
sed -r 's/ /\t/g'  alldepths.tab | sed -r 's/.depth.tab//g' | sort -k1,1n > alldepths.final.tab
