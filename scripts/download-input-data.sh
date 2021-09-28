#!/bin/bash
# Download input data from DDS
# First argument is a folder to download data to
# Second argument is the name of the remote project to download from

DESTINATION=$1
PROJECT=$2

if [ -z "$DESTINATION" ]
then
   echo "ERROR: Missing required destination folder."
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

ddsclient download -p $PROJECT $DESTINATION
