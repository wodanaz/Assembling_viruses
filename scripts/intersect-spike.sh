#!/usr/bin/env bash
# Intersect vcf with spike.bed - Array Job
#
#SBATCH --job-name=ev-intersect-spike-for-depth
#
# Required First Argument: file containing a list of vcf files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

# check positions at the beginning of the line up to the first word boundary in the spike region
POSITIONS_REGEX='^(22812|22813|22917|23012|23063|23403|23592|23593)\b'

root=`basename $VCFFILE .gatk.filt.vcf`
intersectBed -a $VCFFILE -b spike.bed > $EVDIR/$root.spike.vcf

# turn off error checking since grep exits with an error if no items are found
set +e

# create file for the depth compiler
awk '{ print  $2 "\t" $10}' $EVDIR/$root.spike.vcf  | grep -E $POSITIONS_REGEX | sed -r 's/:/\t/g' | awk '{ print  $1 "\t" $3}'> $EVDIR/$root.depth.tab

# create file for the genotype compiler
awk '{ print  $2 "\t" $5}' $EVDIR/$root.spike.vcf  | grep -E $POSITIONS_REGEX > $EVDIR/$root.spike.tab

# make sure they sbatch script doesn't fail if grep finds no items
true
