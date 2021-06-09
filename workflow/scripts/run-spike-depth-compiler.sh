#!/bin/bash
# Compile all tab tables into one for depth
# - tab files to process
# Required environment variable PERLSCRIPT that points to depth compiler

# stop if a command fails (non-zero exit status)
set -e

perl $PERLSCRIPT "$@" > spike_depths.tab
cp spike_depths.tab spike_depths.backup.tab
sed -r 's/ /\t/g' spike_depths.tab | sed -r 's/.depth.tab//g' | sort -k1,1n  > results/spike_depths.final.tab
