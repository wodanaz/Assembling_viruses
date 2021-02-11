#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file - consensus
#
#SBATCH --job-name=ev-bcftools-concensus

# stop if a command fails (non-zero exit status)
set -e

module load bcftools/1.10.2-fasrc01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .gatk.filt.vcf`
bcftools view -Oz $VCFFILE > $VCFFILE.gz
bcftools index $VCFFILE.gz
bcftools norm -f $GENOME $VCFFILE.gz -Ob -o $root.norm.bcf
bcftools consensus -m $root.merged.bed -f $GENOME  -p ${root}_ -s ${root} -H A $VCFFILE.gz > $root.masked.fasta
