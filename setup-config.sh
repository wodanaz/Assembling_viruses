#!/usr/bin/env bash
set -e

echo "This script creates config files allowing single node jobs using scratch storage"

read -p 'Enter target nodename: ' NODENAME
read -p 'Enter scratch base directory: ' SCRATCHDIR
read -p 'Enter account: ' ACCOUNT
read -p 'Enter partition: ' PARTITION

echo ""
echo "Miniconda must be installed at $SCRATCHDIR/miniconda3 with snakemake and dds environments"
echo "The pipeline conda environments will be installed at $SCRATCHDIR/conda"
echo ""
read -p 'Proceed (y/n): ' PROCEED

if [ "$PROCEED" != "y" ]
then
    echo "Exiting without making changes."
    exit 1
fi
PROFILEPATH=$(pwd)/smprofile

echo "Creating config.sh"

cat <<EOT > config.sh
export USE_MODULES=N
export SM_CONDA_PREFIX=$SCRATCHDIR/conda
export SBATCH_FLAGS="--nodelist=$NODENAME"
export SBATCH_PARTITION="$PARTITION"
export SBATCH_ACCOUNT="$ACCOUNT"
export ACTIVATE_CONDA_PATH=$SCRATCHDIR/miniconda3/bin/conda
export SNAKEMAKE_PROFILE=$PROFILEPATH
EOT

echo "Creating snakemake profile config.yaml"
cat <<EOT > smprofile/config.yaml
restart-times: 3
jobscript: "slurm-jobscript.sh"
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --parsable
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=snakejob.{rule}.{jobid}.sh
    --partition=$PARTITION
    --nodelist=$NODENAME
    --account=$ACCOUNT
default-resources:
  - mem_mb=2000
cluster-status: "status-sacct.sh"
max-jobs-per-second: 50
max-status-checks-per-second: 10
local-cores: 1
latency-wait: 60
EOT

echo "Done"
