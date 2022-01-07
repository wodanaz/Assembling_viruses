#!/usr/bin/env bash
set -e

ShowHelp()
{
   # Display Help
   echo "Runs a Slurm pipeline determining escape variants in fastq.gz files, optionally staging data in and out."
   echo
   echo "usage: $0 -d datadir -i inputproject [-g genome] [-j numjobs] [-e email] [-s] [-S] [-k]"
   echo "options:"
   echo "-d datadir       directory used to hold input and output files - required"
   echo "-o outdir        specify where to store output files - optional overrides datadir for output location"
   echo "-i inputproject  project name to download - required"
   echo "-g genome        *.fasta genome to use - defaults to resources/NC_045512.fasta"
   echo "-j numjobs       number of array jobs to run in parallel - defaults to 10"
   echo "-e email         email address to notify on pipeline completion - defaults to empty(no email sent)"
   echo "-s               skip download input data step"
   echo "-S               skip upload output data step"
   echo "-k               keep input/intermediate data directories"
   echo ""
   echo "NOTE: The genome and datadir must be shared across the slurm cluster."
   echo ""
}


###################################
### Read command line arguments ###
###################################


# set default argument values
export GENOME=$(readlink -e resources/NC_045512.fasta)
export OUTDIR=$(pwd)
export DELETE_INT_DIRS=Y
export PROJECTNAME=
export SPIKEBED=$(readlink -e resources/spike.bed)
export NUM_CORES=10
export SNAKEMAKE_PROFILE=$(readlink -e slurm/)
export DOWNLOAD_INPUT_DATA=Y
export UPLOAD_OUTPUT_DATA=Y
export FOREGROUND_MODE=N
export SNAKEMAKE_REPORT=Y

while getopts "g:d:o:i:m:w:j:D:sSke:" OPTION; do
    case $OPTION in
    g)
        export GENOME=$(readlink -e $OPTARG)
        ;;
    d)
        export DATADIR=$(readlink -m $OPTARG)
        ;;
    o)
        OUTPUT_DATADIR=$(readlink -m $OPTARG)
        ;;
    i)
        export PROJECTNAME=$OPTARG
        ;;
    j)
        export NUM_CORES=$OPTARG
        ;;
    e)
        EMAIL=$OPTARG
        ;;
    k)
        export DELETE_INT_DIRS=N
        ;;
    s)
        export DOWNLOAD_INPUT_DATA=N
        ;;
    S)
        export UPLOAD_OUTPUT_DATA=N
        ;;
    esac
done

# check that the config.sh has been setup
if [ -f config.sh ]
then
   source config.sh
else
   echo "ERROR: Missing config.sh config file."
   echo "To fix run:"
   echo "  cp example-config.sh config.sh"
   echo ""
   exit 1
fi

####################################
### Check command line arguments ###
####################################

# set snakemake conda directory if not set from config.sh
if [ -z "$SM_CONDA_PREFIX" ]
then
    export SM_CONDA_PREFIX=$DATADIR/conda
fi

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
if [ -z "$PROJECTNAME" ]
then
   echo "ERROR: Missing required '-i inputproject' argument."
   echo ""
   ShowHelp
   exit 1
fi

# enable USE_MODULES if not explicitly turned off by config.sh
if [ "$USE_MODULES" != "N" ]
then
   export USE_MODULES=Y
fi

# declare directory names
export INPUTDIR=$DATADIR/input/$PROJECTNAME
export SNAKEMAKE_DIR=$DATADIR/output/${PROJECTNAME}_snakemake
if [ -z "$OUTPUT_DATADIR" ]
then
    export OUTDIR="$DATADIR/output/$PROJECTNAME"
else
    export OUTDIR="$OUTPUT_DATADIR/$PROJECTNAME"
fi
export LOGDIR="$OUTDIR/logs"

# make output logs directory if necessary
mkdir -p $LOGDIR

###########################
### Run Escape Variants ###
###########################

SBATCH_FLAGS="$SBATCH_FLAGS --output=$LOGDIR/run-staging-pipeline-%j.out"
if [ ! -z "$EMAIL" ]
then
   echo "Emailing $EMAIL on pipeline completion."
   SBATCH_FLAGS="$SBATCH_FLAGS --mail-type=END --mail-user=$EMAIL"
fi

if [ "$FOREGROUND_MODE" = "Y" ]
then
   bash scripts/run-staging-pipeline.sh | tee $LOGDIR/run-staging-pipeline.out
else
   JOBID=$(sbatch --parsable ${SBATCH_FLAGS} scripts/run-staging-pipeline.sh)
   echo "Submitted batch job $JOBID"
   echo "To monitor main log run:"
   echo "tail -f $LOGDIR/run-staging-pipeline-$JOBID.out"
fi

