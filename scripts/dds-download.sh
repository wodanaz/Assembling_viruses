#!/bin/bash
# Download input data from DDS

module load ddsclient/3.3.0-gcb01
ddsclient download -p $PROJECT $DESTINATION
