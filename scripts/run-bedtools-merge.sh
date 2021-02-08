#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file - Array Job
#
#SBATCH --job-name=ev-concensus-merge
#SBATCH --output=logs/ev-concensus-merge-%j-%a.out

# stop if a command fails (non-zero exit status)
set -e

module load bedtools2/2.25.0-fasrc01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BEDFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BEDFILE .depth.bed`
awk '{ if ( $3 == 0 )  print "MT246667" "\t" $2 "\t" $2 }' ${BEDFILE} > $root.mask.bed
bedtools merge -i $root.mask.bed | awk '{ print $1 "\t" $2 + 1   "\t" $3 - 1 }' > $root.merged.bed
