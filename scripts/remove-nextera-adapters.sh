#!/bin/bash
# Remove Nextera Adapters - Array Job
#
#SBATCH --job-name=ev-rna-ary
#SBATCH --output=logs/ev-rna-%j-%a.out

# stop if a command fails (non-zero exit status)
set -e

module load cutadapt/2.3-gcb01
module load TrimGalore/0.6.5-fasrc01

# Determine the file to process in $FILENAMES_FILE on SLURM_ARRAY_TASK_ID
READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

trim_galore --fastqc --nextera $READFILE
