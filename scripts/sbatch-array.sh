#!/usr/bin/env bash
# Runs sbatch array job.
# Usage: sbatch-array.sh <SCRIPTNAME> <FILENAMES_FILE>
# The SCRIPTNAME should be a sbatch array job script to run.
# The FILENAMES_FILE is a file containing a list of files to process.
# A job will be run for each line in FILENAMES_FILE.
# The sbatch array job should read FILENAMES_FILE based on SLURM_ARRAY_TASK_ID to determine the filename to process.
# eg. READFILE=$(awk NR==$SLURM_ARRAY_TASK_ID $FILENAMES_FILE)

# Controls how many jobs are run at the same time.
MAX_ARRAY_JOBS=4

SCRIPTNAME=$1
# export FILENAMES_FILE for use by array job scripts
export FILENAMES_FILE=$2

NUM_FILES=$(cat ${FILENAMES_FILE} | wc -l)
echo "Processing ${NUM_FILES} files from ${FILENAMES_FILE} using ${SCRIPTNAME}"
# create an array job for each file but only allow MAX_ARRAY_JOBS processes at a time
sbatch --wait --array=1-${NUM_FILES}%${MAX_ARRAY_JOBS} $SCRIPTNAME

