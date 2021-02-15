#!/usr/bin/env bash
set -e

ShowHelp()
{
   # Display Help
   echo "Runs a Slurm pipeline determining escape variants in fastq.gz files, staging data from/to DDS."
   echo
   echo "usage: $0 -g genome -d datadir -i inputproject"
   echo "options:"
   echo "-g genome        *.fasta genome to use - required"
   echo "-d datadir       directory used to hold input and output files - required"
   echo "-i inputproject  project name to download - required"
   echo "-s               runs surveillance mode - default is run experimental mode"
   echo ""
   echo "NOTE: The input genome must first be indexed by running ./setup-variants-pipeline.sh."
   echo "NOTE: The genome and datadir must be shared across the slurm cluster."
   echo ""
}

REV_ARGS=""
while getopts "g:d:i:s" OPTION; do
    case $OPTION in
    g)
        export GENOME=$(readlink -e $OPTARG)
        ;;
    d)
        export DATADIR=$(readlink -e $OPTARG)
        ;;
    i)
        export DDS_INPUT_PROJECT=$OPTARG
        ;;
    s)
        REV_ARGS="$REV_ARGS -s"
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
if [ -z "DATADIR" ]
then
   echo "ERROR: Missing required '-d datadir' argument."
   echo ""
   ShowHelp
   exit 1
fi
if [ -z "DDS_INPUT_PROJECT" ]
then
   echo "ERROR: Missing required '-i inputproject' argument."
   echo ""
   ShowHelp
   exit 1
fi


# declare directory names
INPUT_DATADIR=$DATADIR/input
INPUT_PROJECT_DIR=$INPUT_DATADIR/$DDS_INPUT_PROJECT
OUTPUT_DATADIR=$DATADIR/output
OUTPUT_RESULTS_DIR=$OUTPUT_DATADIR/$DDS_INPUT_PROJECT


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
bash scripts/dds-download.sh
echo ""


echo "Running escape variants pipeline - logs at $OUTPUT_RESULTS_DIR/logs"
# run pipeline
SBATCH_FLAGS="--wait" ./run-escape-variants.sh \
  $REV_ARGS \
  -g $GENOME \
  -i $INPUT_PROJECT_DIR \
  -o $OUTPUT_RESULTS_DIR \
  -l $OUTPUT_RESULTS_DIR/logs
echo ""


# upload results to DDS
export PROJECT=${DDS_INPUT_PROJECT}_results
export SOURCE="$OUTPUT_RESULTS_DIR/*"
echo "Uploading $DESTINATION to $PROJECT"
bash scripts/dds-upload.sh
echo ""


echo "Deleting $INPUT_PROJECT_DIR"
rm -rf $INPUT_PROJECT_DIR


echo "DDS Escape Variants Done"
