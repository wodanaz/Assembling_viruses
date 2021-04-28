#!/bin/bash
# Add sme info for the read groups
#
#SBATCH --job-name=ev-update-rg
#SBATCH --mem 4000
#
# Required First Argument: file containing a list of *dedup.bam files to process
# Required Environment Variables:
#  EVDIR - working directory
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .dedup.bam`;

picard -Xmx7g AddOrReplaceReadGroups -Djava.io.tmpdir=/data/covid19lab/tmp I=$BAMFILE O=$EVDIR/$root.bam2 RGSM=$root RGPU=unit1 RGLB=lib_${root} RGPL=ILLUMINA
