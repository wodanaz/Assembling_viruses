#!/bin/bash
# Compile all tab tables into one for depth
#
#SBATCH --job-name=ev-run-spike-depth-compiler
#
# Required Environment Variables:
#  EVDIR - working directory containing *.depth.tab files

# stop if a command fails (non-zero exit status)
set -e

perl depth_compiler.pl $EVDIR/*.depth.tab > $EVDIR/spike_depths.tab
cp $EVDIR/spike_depths.tab $EVDIR/spike_depths.backup.tab
sed -r 's/ /\t/g' $EVDIR/spike_depths.tab | sed -r 's/.depth.tab//g' | sort -k1,1n  > $EVDIR/spike_depths.final.tab
