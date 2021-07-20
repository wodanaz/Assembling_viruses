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
echo ""

if [ "$DOWNLOAD_INPUT_DATA" == "Y" ]
then
    # download project from DDS
    echo "Downloading input data to $INPUTDIR"
    srun --output="${LOGDIR}/download-input-data-%j.out" scripts/download-input-data.sh $INPUTDIR $PROJECTNAME
    echo ""
fi

echo "Running escape variants pipeline - logs at $LOGDIR"
srun --output="${LOGDIR}/escape-variants-pipeline-%j.out" scripts/escape-variants-pipeline.sh
echo ""

if [ "$UPLOAD_OUTPUT_DATA" == "Y" ]
then
    # upload results to DDS
    echo "Uploading results from $OUTDIR"
    srun --output="${LOGDIR}/upload-output-data-%j.out" scripts/upload-output-data.sh $OUTDIR ${PROJECTNAME}_results
    echo ""
fi

echo "Deleting $INPUTDIR"
rm -rf $INPUTDIR

echo "Escape Variants Done"
