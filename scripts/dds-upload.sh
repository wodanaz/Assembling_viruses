#!/bin/bash
# Upload results to DDS

conda activate dds
ddsclient upload -p $PROJECT $SOURCE
