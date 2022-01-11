#!/usr/bin/env bash
# Script that will optionally stage data and run the escape variants pipeline.
# Required environment variables
# - PROJECTNAME - Name of the project to download from (uploads to "${PROJECTNAME}_results")
# - DOWNLOAD_INPUT_DATA - when "Y" download input data into $INPUTDIR
# - UPLOAD_OUTPUT_DATA - when "Y" upload results into $INPUTDIR
# - INPUTDIR - directory to hold input data
# - OUTDIR - directory to hold output data
# - LOGDIR - directory to hold log files
#
#SBATCH --job-name=ev-staging-pipeline

# stop if a command fails (non-zero exit status)
set -e

echo "Escape Variants Starting"
date
echo ""

# make input results directory if necessary
mkdir -p INPUTDIR
# make output results directory if necessary
mkdir -p $OUTDIR
# create temp snakemake working directory <datadir>/output/<projectname>_snakemake
mkdir -p $SNAKEMAKE_DIR
# create directory to hold conda environments if necessary
mkdir -p $SM_CONDA_PREFIX

# Activate a conda shell if configured by config.sh
if [ ! -z "$ACTIVATE_CONDA_PATH" ]
then
   echo "Activate $ACTIVATE_CONDA_PATH."
   eval "$($ACTIVATE_CONDA_PATH shell.bash hook)"
fi

if [ "$DOWNLOAD_INPUT_DATA" == "Y" ]
then
    # download project from DDS
    echo "Downloading input data to $INPUTDIR"
    if [ "$FOREGROUND_MODE" = "Y" ]
    then
        bash scripts/download-input-data.sh $INPUTDIR $PROJECTNAME
    else
        srun --output="${LOGDIR}/download-input-data-%j.out" scripts/download-input-data.sh $INPUTDIR $PROJECTNAME
    fi
    echo ""
fi

echo "Running escape variants pipeline"
date
./scripts/escape-variants-pipeline.sh 2>&1
echo ""

echo "Running escape variants pipeline - DONE"
date
echo ""

if [ "$UPLOAD_OUTPUT_DATA" == "Y" ]
then
    # upload results to DDS
    echo "Uploading results from $OUTDIR"
    if [ "$FOREGROUND_MODE" = "Y" ]
    then
        bash scripts/upload-output-data.sh $OUTDIR ${PROJECTNAME}_results
    else
        srun --output="${LOGDIR}/upload-output-data-%j.out" scripts/upload-output-data.sh $OUTDIR ${PROJECTNAME}_results
    fi
    echo ""
fi

if [ "$DELETE_INT_DIRS" == "Y" ]
    then
    echo "Deleting temporary $SNAKEMAKE_DIR directory"
    rm -rf $SNAKEMAKE_DIR

    echo "Deleting $INPUTDIR"
    rm -rf $INPUTDIR
fi


echo "Escape Variants Done"
