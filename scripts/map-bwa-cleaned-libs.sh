#!/bin/bash
# Map using BWA with the cleaned libraries - Array Job
#
#SBATCH --job-name=ev-map-bwa-ary
#
# Required First Argument: file containing a list of *_trimmed.fq.gz files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

module load bwa/0.7.12-gcb01

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=$(basename $READFILE _trimmed.fq.gz)
root2=$(basename $root _R1_001)
root3=$(echo ${root} | cut -d'_' -f 1)

bwa mem $GENOME -R "@RG\tID:ID_${root3}\tPU:PU_${root3}\tSM:${root3}\tLB:${root}" $READFILE > $EVDIR/$root3.sam
