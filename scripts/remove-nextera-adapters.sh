#!/bin/bash
# Remove Nextera Adapters
# Required positional arguments:
# - file to trim in fastq.gz format

# stop if a command fails (non-zero exit status)
set -e

READFILE=$1

trim_galore --fastqc --nextera $READFILE
