#!/usr/bin/env bash
# Intersect vcf with a spike bed file
# - vcf file to process
# - spike file in bed format

# stop if a command fails (non-zero exit status)
set -e

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$1
SPIKEFILE=$2

# check positions at the beginning of the line up to the first word boundary in the spike region
POSITIONS_REGEX='^(22812|22813|22917|23012|23063|23403|23592|23593)\b'

root=`basename $VCFFILE .gatk.filt.vcf`
intersectBed -a $VCFFILE -b $SPIKEFILE > $root.spike.vcf

# turn off error checking since grep exits with an error if no items are found
set +e

# create file for the depth compiler
awk '{ print  $2 "\t" $10}' $root.spike.vcf  | grep -E $POSITIONS_REGEX | sed -r 's/:/\t/g' | awk '{ print  $1 "\t" $3}'> $root.depth.tab

# create file for the genotype compiler
awk '{ print  $2 "\t" $5}' $root.spike.vcf  | grep -E $POSITIONS_REGEX > $root.spike.tab
