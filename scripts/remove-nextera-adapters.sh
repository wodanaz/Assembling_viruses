#!/bin/bash
# Remove Nextera Adapters
# Required positional arguments:
# - file to trim in fastq.gz format

# stop if a command fails (non-zero exit status)
set -e

READFILE=$1
OUTFILE=$2

trim_galore --fastqc --nextera $READFILE

BASEREAD=$(basename $READFILE)
mv "${BASEREAD/.fastq.gz/_trimmed.fq.gz}" $OUTFILE
