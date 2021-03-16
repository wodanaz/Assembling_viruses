#!/bin/bash
# Make final consensus fasta sequence using the SNPs in the vcf file and trims first and last 25 bp - consensus
#
#SBATCH --job-name=ev-bcftools-concensus
#
# Required First Argument: file containing a list of *gatk.filt.vcf files to process
# Required Environment Variables:
#  EVDIR - working directory
#  FILENAMES_FILE - list of *gatk.filt.vcf files to process
#  SLURM_ARRAY_TASK_ID - number specifying which file in $1 (FILENAMES_FILE) to process

# stop if a command fails (non-zero exit status)
set -e

# The first argument is a file containing filenames to process
FILENAMES_FILE=$1
# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
VCFFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $VCFFILE .gatk.filt.vcf`
bcftools view -Oz $VCFFILE > $VCFFILE.gz
bcftools index $VCFFILE.gz
bcftools norm -f $GENOME $VCFFILE.gz -Ob -o $EVDIR/$root.norm.bcf
bcftools consensus -m $EVDIR/$root.merged.bed -f $GENOME  -p ${root}_ -s ${root} -H A $VCFFILE.gz > $EVDIR/$root.masked.fasta
seqtk trimfq -b 25 -e 25  $EVDIR/$root.masked.fasta > $EVDIR/$root.cleaned.fasta
