#!/bin/bash
# BAM TO VCF USING BCFTOOLS TO CREATE A PRELIMINAR VCF FILE
#
#SBATCH --job-name=ev-gatk-bam-to-vcf
#SBATCH --output=logs/ev-gbam-to-vcf-%j-%a.out
#SBATCH --mem 10G

# stop if a command fails (non-zero exit status)
set -e

module load samtools/1.10-gcb01
module load bcftools/1.10.2-fasrc01

# Determine the file to process in $FILENAMES_FILE based on SLURM_ARRAY_TASK_ID
BAMFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

root=`basename $BAMFILE .dedup.bam`;

bcftools mpileup -Ou -f sars_cov_2.fasta $BAMFILE --annotate FORMAT/DPR > $root.bcf
bcftools call -vm --ploidy 1 $root.bcf > $root.raw.vcf
