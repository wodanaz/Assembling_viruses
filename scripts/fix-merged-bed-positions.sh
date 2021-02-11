#!/usr/bin/env bash
# Fix merged bed positions
# For some ridiculous reason, bedtools messes up the first and last base when merging and they must be edited
# Here, a general rule to change the positions affected
#
#SBATCH --job-name=ev-fixbed

# stop if a command fails (non-zero exit status)
set -e

sed -ri 's/MT246667\t29872\t29870/MT246667\t29870\t29871/g' *merged.bed
sed -ri 's/MT246667\t2\t0/MT246667\t0\t1/g' *merged.bed
sed -ri 's/MT246667\t1\t/MT246667\t0\t/g' *merged.bed
