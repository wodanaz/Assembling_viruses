#!/usr/bin/env bash
set -e

ShowHelp()
{
   # Display Help
   echo "Runs a Slurm pipeline determining escape variants in fastq.gz files, staging data from/to DDS."
   echo
   echo "usage: $0 -g genome -d datadir -i inputproject -m mode [-w workdir] [-j numjobs] [-D datetab] [-k]"
   echo "options:"
   echo "-g genome        *.fasta genome to use - required"
   echo "-d datadir       directory used to hold input and output files - required"
   echo "-i inputproject  project name to download - required"
   echo "-m mode          run mode: 'h' (hospital surveillance) or 'c' (campus surveillance ) or 'e' (experimental) - required"
   echo "-w workdir       directory that will hold a tempdir - defaults to current directory"
   echo "-j numjobs       number of array jobs to run in parallel - defaults to 4"
   echo "-D datetab       date.tab file used to create supermetadata.tab - defaults to skipping the supermetadata step"
   echo "-k               Keep working directory. When passed the working directory will not be deleted"
   echo ""
   echo "NOTE: The genome and datadir must be shared across the slurm cluster."
   echo ""
}

export NUM_CORES=10
export SNAKEMAKE_PROFILE=$(readlink -e slurm/)
DATETAB_OPT=""
GENOME=$(readlink -e resources/NC_045512.fasta)
SPIKEBED=$(readlink -e resources/spike.bed)
CLEANUP_WORKDIR=Y
while getopts "g:d:i:m:w:j:D:k" OPTION; do
    case $OPTION in
    g)
        GENOME=$(readlink -e $OPTARG)
        ;;
    d)
        DATADIR=$(readlink -f $OPTARG)
        ;;
    m)
        RUN_MODE=$OPTARG
        ;;
    i)
        DDS_INPUT_PROJECT=$OPTARG
        ;;
    j)
        export NUM_CORES=$OPTARG
        ;;
    w)
        WORKDIR=$(readlink -e $OPTARG)
        ;;
    D)
        DATETAB_OPT="datetab: $(readlink -e $OPTARG)"
        ;;
    k)
        CLEANUP_WORKDIR=N
        ;;
    esac
done


# check required arguments
if [ -z "$GENOME" ]
then
   echo "ERROR: Missing required '-g genome' argument."
   echo ""
   ShowHelp
   exit 1
fi
if [ -z "$DATADIR" ]
then
   echo "ERROR: Missing required '-d datadir' argument."
   echo ""
   ShowHelp
   exit 1
fi
if [ -z "$DDS_INPUT_PROJECT" ]
then
   echo "ERROR: Missing required '-i inputproject' argument."
   echo ""
   ShowHelp
   exit 1
fi

case $RUN_MODE in
  c)
    RUN_MODE="campus"
    ;;
  h)
    RUN_MODE="hospital"
    ;;
  e)
    RUN_MODE="experimental"
    ;;
  *)
    echo "ERROR: Invalid '-m mode' argument: $RUN_MODE"
    exit 1
    ;;
esac


# declare directory names
INPUT_DATADIR=$DATADIR/input
INPUT_PROJECT_DIR=$INPUT_DATADIR/$DDS_INPUT_PROJECT/
export WORKDIR=$DATADIR/work/$DDS_INPUT_PROJECT/
OUTPUT_DATADIR=$DATADIR/output
OUTPUT_RESULTS_DIR=$OUTPUT_DATADIR/$DDS_INPUT_PROJECT
export SM_CONDA_PREFIX=$DATADIR/conda

echo "DDS Escape Variants Starting"
echo ""


# make base input directory if necessary
mkdir -p $INPUT_DATADIR
# make base output directory if necessary
mkdir -p $OUTPUT_DATADIR
# make output results directory if necessary so logs parent directory exists
mkdir -p $OUTPUT_RESULTS_DIR


# download project from DDS
export PROJECT=$DDS_INPUT_PROJECT
export DESTINATION=$INPUT_PROJECT_DIR
echo "Downloading $PROJECT to $DESTINATION"
./scripts/dds-download.sh
echo ""

# create config/config.yaml
mkdir -p $WORKDIR/config
cat <<EOF > $WORKDIR/config/config.yaml
project: $DDS_INPUT_PROJECT
genome: $GENOME
spike: $SPIKEBED
inputdir: $INPUT_PROJECT_DIR
mode: $RUN_MODE
$DATETAB_OPT
EOF
echo "Created snakemake $WORKDIR/config/config.yaml file:"
cat $WORKDIR/config/config.yaml
echo ""

echo "Running escape variants pipeline - logs at $WORKDIR/logs"
./scripts/run-snakemake.sh

echo "Copying results to $OUTPUT_RESULTS_DIR"
cp -r $WORKDIR/results/* $OUTPUT_RESULTS_DIR/.
cp -r $WORKDIR/logs $OUTPUT_RESULTS_DIR/logs

## upload results to DDS
export PROJECT=${DDS_INPUT_PROJECT}_results
export SOURCE="$OUTPUT_RESULTS_DIR/*"
echo "Uploading $DESTINATION to $PROJECT"
bash scripts/dds-upload.sh
echo ""

if [ "$CLEANUP_WORKDIR" == "Y" ]
then
    echo "Deleting $WORKDIR"
    rm -rf $WORKDIR
fi

echo "Deleting $INPUT_PROJECT_DIR"
rm -rf $INPUT_PROJECT_DIR

echo "DDS Escape Variants Done"
