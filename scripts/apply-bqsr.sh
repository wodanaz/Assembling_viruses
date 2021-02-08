#!/bin/bash
# APPLY BQSR (Apply a linear base quality recalibration model trained with the BaseRecalibrator tool)
#
#SBATCH --job-name=ev-apply-bqsr
#SBATCH --output=logs/ev-apply-bqsr-%j-%a.out
#SBATCH --mem 2000

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bam2`
gatk --java-options -Xmx8G  ApplyBQSR -I $BAMFILE -R MT246667.fasta --bqsr-recal-file $root.table  -O $root.bqsr.bam
