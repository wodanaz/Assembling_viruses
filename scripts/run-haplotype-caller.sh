#!/bin/bash
# HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes
#
#SBATCH --job-name=ev-haplotype-caller
#SBATCH --mem 2000
#
# Required First Argument: file containing a list of *bqsr.bam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bqsr.bam`
gatk --java-options "-Djava.io.tmpdir=/data/covid19lab/tmp" HaplotypeCaller -I $BAMFILE -R $GENOME -ploidy 1 -O $EVDIR/$root.gatk.vcf
