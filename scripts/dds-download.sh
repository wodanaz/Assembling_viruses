#!/bin/bash
# Download input data from DDS

module load ddsclient/3.2.0-gcb01
ddsclient download -p $PROJECT $DESTINATION
