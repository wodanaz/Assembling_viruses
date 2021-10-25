#!/bin/bash
# Map using BWA with the cleaned libraries
# Required positional arguments:
# - file in fastq.gz format used to create a sam file
# - indexed genome file used by bwa

# stop if a command fails (non-zero exit status)
set -e

READFILE=$1
GENOME=$2

root=$(basename $READFILE _trimmed.fq.gz)
root2=$(basename $root _R1_001)
root3=$(echo ${root} | cut -d'_' -f 1)

bwa mem -R "@RG\tID:ID_${root3}\tPU:PU_${root3}\tSM:${root3}\tLB:${root}" $GENOME $READFILE > $root3.sam
