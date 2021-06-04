#!/usr/bin/env bash
# Index Genome Reference
# Required positional argument: genome file to index

# stop if a command fails (non-zero exit status)
set -e

GENOME=$1

bwa index -a is $GENOME
samtools faidx $GENOME
