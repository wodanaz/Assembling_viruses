#!/bin/bash
# Add sme info for the read groups
#
#SBATCH --job-name=ev-update-rg
#SBATCH --output=logs/ev-update-rg-%j-%a.out
#SBATCH --mem 1000

# stop if a command fails (non-zero exit status)
set -e

module load picard-tools/2.4.1-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .dedup.bam`;

java -Xmx7g -jar $PICARD_TOOLS_HOME/picard.jar AddOrReplaceReadGroups I=$BAMFILE O=$root.bam2 RGSM=$root RGPU=unit1 RGLB=lib_${root} RGPL=ILLUMINA
