#!/bin/bash
# Upload results to DDS
# First argument is a folder to upload data from
# Second argument is the name of the remote project to upload into

SOURCE=$1
PROJECT=$2

if [ -z "$SOURCE" ]
then
   echo "ERROR: Missing required source folder."
   echo ""
   exit 1
fi

if [ "$USE_MODULES" == "Y" ]
then
    module load ddsclient/3.2.0-gcb01
else
    # make sure conda is setup within this script
    source $CONDA_PREFIX/etc/profile.d/conda.sh
    conda activate dds
fi

ddsclient upload -p $PROJECT $SOURCE/*
