#!/bin/bash
# Base recalibration. First pass of the Base Quality Score Recalibration (BQSR)
# Generates recalibration table based on various user-specified covariates
# (such as read group, reported quality score, machine cycle, and nucleotide context).
#
#SBATCH --job-name=ev-baserecal
#SBATCH --mem 2000

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01
module load tabix/0.2.6-fasrc01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bam2`

tabix -p vcf $root.filt.vcf.gz -f
gatk --java-options -Xmx8G BaseRecalibrator -I $BAMFILE -R $GENOME --known-sites $root.filt.vcf.gz -O $root.table
