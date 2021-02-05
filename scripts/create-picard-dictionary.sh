#!/usr/bin/env bash
# Create a dictionary file for using picard tools
#
#SBATCH --job-name=ev-dict
#SBATCH --output=logs/ev-dict-%j.out

# stop if a command fails (non-zero exit status)
set -e

module load picard-tools/2.4.1-gcb01
module load samtools/1.10-gcb01

java -jar /nfs/software/helmod/apps/Core/picard-tools/2.4.1-gcb01/picard.jar  CreateSequenceDictionary R=MT246667.fasta O=MT246667.dict 
samtools faidx MT246667.fasta
