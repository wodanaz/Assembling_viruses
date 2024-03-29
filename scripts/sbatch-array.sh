#!/usr/bin/env bash
# Runs sbatch array job.
# Usage: sbatch-array.sh <SCRIPTNAME> <FILENAMES_FILE>
# The SCRIPTNAME should be a sbatch array job script to run.
# The FILENAMES_FILE is a file containing a list of files to process.
# A job will be run for each line in FILENAMES_FILE.
# The sbatch array job should read FILENAMES_FILE based on SLURM_ARRAY_TASK_ID to determine the filename to process.
# eg. READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

# MAX_ARRAY_JOBS controls how many jobs are run at the same time.
if [ -z "$MAX_ARRAY_JOBS" ]
then
  MAX_ARRAY_JOBS=4
fi

SCRIPTNAME=$1
FILENAMES_FILE=$2

# create a logfile based on the name script we are running
SCRIPTBASE=$(basename ${SCRIPTNAME} .sh)
SBATCH_OUTPUT_FLAG="--output=${LOGDIR}/${SCRIPTBASE}-%A_%a.out"

NUM_FILES=$(cat ${FILENAMES_FILE} | wc -l)
echo "Processing ${NUM_FILES} file(s) using:"
echo "    ${SCRIPTNAME} ${FILENAMES_FILE}"
# create an array job for each file but only allow MAX_ARRAY_JOBS processes at a time
sbatch --wait --array=1-${NUM_FILES}%${MAX_ARRAY_JOBS} $SBATCH_OUTPUT_FLAG $SCRIPTNAME $FILENAMES_FILE
