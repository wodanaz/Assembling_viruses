#!/bin/bash
# Creates consensus fasta file, runs pangolin, and checks lineage report for specific values 
#
#SBATCH --job-name=ev-pangolin
#
# Required Environment Variables:
#  EVDIR - working directory containing *.fasta files
#  PROJECTNAME - controls output filenames

# stop if a command fails (non-zero exit status)
set -e

# concatenates all new consensus sequences and runs pangolin
cat $EVDIR/*.cleaned.fasta > $EVDIR/consensus_sequences.fasta
pangolin $EVDIR/consensus_sequences.fasta  --outdir $EVDIR --outfile $PROJECTNAME.csv

# stop checking exit status since grep fails when not found
set +e

# search for variants of concern
grep -E 'B\.1\.351|B\.1\.1\.7|P\.1|P\.2|B\.1\.427|B\.1\.429|B\.1\.526' $EVDIR/$PROJECTNAME.csv > $EVDIR/${PROJECTNAME}_lineages_of_concern.csv

# exit with status 0 so slurm will ignoring the grep exit status
exit 0
