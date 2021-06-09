#!/usr/bin/env bash
# Run escape variants snakemake pipeline

# stop if a command fails (non-zero exit status)
set -e

module load Anaconda3/2019.10-gcb02
conda activate snakemake

echo "Running: snakemake --cores $NUM_CORES --directory $WORKDIR --use-conda --conda-prefix $SM_CONDA_PREFIX --profile $SNAKEMAKE_PROFILE"

snakemake \
    --cores $NUM_CORES \
    --directory $WORKDIR \
    --use-conda --conda-prefix $CONDA_PREFIX \
    --profile $SNAKEMAKE_PROFILE
