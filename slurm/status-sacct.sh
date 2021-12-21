#!/usr/bin/env bash

# Check status of Slurm job

jobid="$1"

if [[ "$jobid" == Submitted ]]
then
  echo status-sacct: Invalid job ID: "$jobid" >&2
  echo status-sacct: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk '{print $1}'`

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]
then
  echo running
else
  echo failed
fi
