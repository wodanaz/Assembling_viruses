#!/bin/bash
# 'Deduplicate' or mark PCR duplicates
#
#SBATCH --job-name=ev-markdup
#SBATCH --mem 10G

# stop if a command fails (non-zero exit status)
set -e

module load picard-tools/2.4.1-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bam`

java -Xmx7g -jar $PICARD_TOOLS_HOME/picard.jar MarkDuplicates I=$BAMFILE O=$root.dedup.bam M=$root.metric.txt
