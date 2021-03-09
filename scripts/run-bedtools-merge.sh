#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file - Array Job
#
#SBATCH --job-name=ev-concensus-merge
#
# Required First Argument: file containing a list of *depth.bed files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

module load bedtools2/2.25.0-fasrc01

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BEDFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BEDFILE .depth.bed`

awk '{ if ( $3 == 0 )  print $1 "\t" $2 "\t" $2 + 1 }' ${BEDFILE} | bedtools merge |  awk '{ print $1 "\t" $2 "\t" $3 - 1  }' > $EVDIR/$root.merged.bed
