#!/bin/bash
# HAPLOTYPE CALLER: Call germline SNPs and indels via local re-assembly of haplotypes
#
#SBATCH --job-name=ev-haplotype-caller
#SBATCH --mem 2000

# stop if a command fails (non-zero exit status)
set -e

module load GATK/4.1.3.0-gcb01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .bqsr.bam`
gatk --java-options -Xmx8G HaplotypeCaller -I $BAMFILE -R $GENOME -ploidy 1 -O $root.gatk.vcf