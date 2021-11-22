#!/bin/bash
# print a table with coverage information
# Required positional arguments:
# - files in tab format to process

# stop if a command fails (non-zero exit status)
set -e

# merge input files
cat "$@" > coverage.tab

sort -k1,1 -V coverage.tab | sed 1i"sampleID\treference\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
