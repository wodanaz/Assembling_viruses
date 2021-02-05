#!/usr/bin/env bash
# Index Genome Reference
#
#SBATCH --job-name=ev-igv
#SBATCH --output=logs/ev-igv-%j.out

# stop if a command fails (non-zero exit status)
set -e

module load bwa/0.7.12-gcb01

bwa index -a is MT246667.fasta

