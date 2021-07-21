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

module load ddsclient/3.2.0-gcb01
ddsclient upload -p $PROJECT $SOURCE/*
