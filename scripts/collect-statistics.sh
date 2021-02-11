#!/bin/bash
# Collect statistics: Produces a summary of alignment metrics from a SAM or BAM file
#
#SBATCH --job-name=ev-collect-stats
#SBATCH --mem 2000

# stop if a command fails (non-zero exit status)
set -e

module load picard-tools/2.4.1-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bqsr.bam`;

java -Xmx7g -jar $PICARD_TOOLS_HOME/picard.jar CollectAlignmentSummaryMetrics R=$GENOME I=$BAMFILE O=$root.stat.txt
