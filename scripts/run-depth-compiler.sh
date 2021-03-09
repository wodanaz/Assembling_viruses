#!/bin/bash
# Compile all tab tables into one for depth
#
#SBATCH --job-name=ev-depth-compiler
#
# Required Environment Variables:
#  EVDIR - working directory containing *.depth.tab files

# stop if a command fails (non-zero exit status)
set -e

module load perl/5.10.1-fasrc04

perl depth_compiler.pl $EVDIR/*.depth.tab > $EVDIR/alldepths.tab
cp $EVDIR/alldepths.tab $EVDIR/alldepths.backup.tab
sed -r 's/ /\t/g'  $EVDIR/alldepths.tab | sed -r 's/.depth.tab//g' | sort -k1,1n > $EVDIR/alldepths.final.tab
