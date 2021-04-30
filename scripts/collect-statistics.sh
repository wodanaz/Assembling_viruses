#!/bin/bash
# Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file
#
#SBATCH --job-name=ev-collect-stats
#SBATCH --mem 4000
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

root=`basename $BAMFILE .bqsr.bam`;

picard -Xmx7g -Djava.io.tmpdir=$EVDIR CollectAlignmentSummaryMetrics R=$GENOME I=$BAMFILE O=$EVDIR/$root.stat.txt
