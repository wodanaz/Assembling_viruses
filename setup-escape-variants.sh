#!/usr/bin/env bash
# Performs one time genome setup required by run-escape-variants.sh
#
# Expects a MT246667.fasta file to be present
set -e

# create logs directory if necessary
mkdir -p logs


echo "Setup for escape-variants-pipeline - Starting"


echo "Setup part 1 - Index Genome Reference"
sbatch --wait scripts/index-reference-genome.sh
echo "Setup part 1 - Done"
echo ""


echo "Setup part 2 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-picard-dictionary.sh
echo "Setup part 2 - Done"
echo ""


echo "Setup for escape-variants-pipeline - Done"
