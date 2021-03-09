#!/bin/bash
# Upload results to DDS

module load ddsclient/3.2.0-gcb01
ddsclient upload -p $PROJECT $SOURCE
