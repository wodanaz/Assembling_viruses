#!/bin/bash
# Creates consensus fasta file, runs pangolin, and checks lineage report for specific values 
# - files in fasta format to process
# Required environment variable PROJECTNAME that points to depth compiler

# stop if a command fails (non-zero exit status)
set -e

# fetch pangolin details
pangolin --update

# concatenates all new consensus sequences and runs pangolin
cat "$@" > consensus_sequences.fasta
pangolin consensus_sequences.fasta  --outdir results --outfile $PROJECTNAME.csv
