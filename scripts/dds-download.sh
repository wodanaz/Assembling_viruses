#!/bin/bash
# Download input data from DDS

conda activate dds
ddsclient download -p $PROJECT $DESTINATION
