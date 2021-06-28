#!/usr/bin/env bash
# Performs one time genome setup required by run-escape-variants.sh
#
set -e

while getopts "g:" OPTION; do
    case $OPTION in
    g)
        export GENOME=$OPTARG
        ;;
    esac
done

if [ -z "$GENOME" ]
then
   echo "ERROR: You must specify a genome via the '-g <GENOME.fasta>' argument."
   echo "Example: $0 -g MT246667.fasta"
   exit 1
fi

# create logs directory if necessary
mkdir -p logs


echo "Setup for escape-variants-pipeline - Starting"

echo "Activating 'escapevariants' conda environment"

# set default value for ANACONDAMODULE to "Anaconda3/2019.10-gcb02"
ANACONDAMODULE="${ANACONDAMODULE-Anaconda3/2019.10-gcb02}"

conda activate escapevariants
pangolin --update

echo ""

echo "Setup part 1 - Index Genome Reference"
sbatch --wait scripts/index-reference-genome.sh
echo "Setup part 1 - Done"
echo ""


echo "Setup part 2 - Create a dictionary file for using picard tools"
sbatch --wait scripts/create-picard-dictionary.sh
echo "Setup part 2 - Done"
echo ""


echo "Setup for escape-variants-pipeline - Done"
