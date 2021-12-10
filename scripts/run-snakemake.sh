#!/usr/bin/env bash
# Run escape variants snakemake pipeline

# stop if a command fails (non-zero exit status)
set -e

if [ "$USE_MODULES" == "Y" ]
then
  # set default value for ANACONDAMODULE to "Anaconda3/2019.10-gcb02"
  ANACONDAMODULE="${ANACONDAMODULE-Anaconda3/2019.10-gcb02}"
  module load $ANACONDAMODULE
else
  # make sure conda is setup within this script
  source $CONDA_PREFIX/etc/profile.d/conda.sh
fi

conda activate snakemake

echo "Running: snakemake --cores $NUM_CORES --directory $SNAKEMAKE_DIR --use-conda --conda-prefix $SM_CONDA_PREFIX --profile $SNAKEMAKE_PROFILE"

snakemake \
    --cores $NUM_CORES \
    --directory $SNAKEMAKE_DIR \
    --use-conda --conda-prefix $SM_CONDA_PREFIX \
    --profile $SNAKEMAKE_PROFILE

if [ "$SNAKEMAKE_REPORT" == "Y" ]
then
    REPORT_NAME=results/report-$(date +"%Y-%m-%d-%T").html
    echo "Creating snakemake report $REPORT_NAME"
    snakemake --report $REPORT_NAME --directory $SNAKEMAKE_DIR
fi
