#!/usr/bin/env bash
# Investigating escape mutations in SARS-CoV-2
#

# stop if a command fails (non-zero exit status)
set -e

# show an error message when a step fails
trap "echo 'ERROR: Script failed'" ERR

# record the hash of the latest commit of this script
export EV_GIT_COMMIT=$(git log --pretty="%H" -n 1)

echo "Pipeline - Starting - Version $EV_GIT_COMMIT"
echo ""

echo "Current directory is $(pwd)"
echo ""

echo "Running using snakemake directory $SNAKEMAKE_DIR"
echo ""

# create config/config.yaml
mkdir -p $SNAKEMAKE_DIR/config
cat <<EOF > $SNAKEMAKE_DIR/config/config.yaml
project: $PROJECTNAME
genome: $GENOME
spike: $SPIKEBED
inputdir: $INPUTDIR
readsuffix: $READ_SUFFIX
EOF
echo "Created snakemake $SNAKEMAKE_DIR/config/config.yaml file:"
cat $SNAKEMAKE_DIR/config/config.yaml
echo ""

echo "Running escape variants pipeline - logs at $SNAKEMAKE_DIR/logs"
./scripts/run-snakemake.sh

echo "Copying results to $OUTDIR"
cp -r $SNAKEMAKE_DIR/results/* $OUTDIR/.
tar -C $SNAKEMAKE_DIR/logs -czf $OUTDIR/logs.tar.gz .

echo ""
echo "Pipeline - Done - Results in $OUTDIR"
