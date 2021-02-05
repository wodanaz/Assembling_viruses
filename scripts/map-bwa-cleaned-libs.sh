#!/bin/bash
# Map using BWA with the cleaned libraries - Array Job
#
#SBATCH --job-name=ev-map-ary
#SBATCH --output=logs/ev-map-%j-%a.out

# stop if a command fails (non-zero exit status)
set -e

module load bwa/0.7.12-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=$(basename $READFILE _trimmed.fq.gz)
root2=$(basename $root _R1_001)
root3=$(echo ${root} | cut -d'_' -f 1)

bwa mem MT246667.fasta -R "@RG\tID:ID_${root3}\tPU:PU_${root3}\tSM:${root3}\tLB:${root}" $READFILE > $root3.sam

