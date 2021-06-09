#!/bin/bash
# Creates consensus fasta file, runs pangolin, and checks lineage report for specific values 
# - files in fasta format to process
# Required environment variable PROJECTNAME that points to depth compiler

# stop if a command fails (non-zero exit status)
set -e
set -x

# fetch pangolin details
pangolin --update

# concatenates all new consensus sequences and runs pangolin
cat "$@" > consensus_sequences.fasta
pangolin consensus_sequences.fasta  --outdir results --outfile $PROJECTNAME.csv

# stop checking exit status since grep fails when not found
set +e

# search for variants of concern
grep -E 'B\.1\.351|B\.1\.1\.7|P\.1|P\.2|B\.1\.427|B\.1\.429|B\.1\.526' $PROJECTNAME.csv > results/${PROJECTNAME}_lineages_of_concern.csv

# exit with status 0 so slurm will ignoring the grep exit status
exit 0
