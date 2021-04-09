#!/usr/bin/env bash
# Escape Variants
#
# Runs escape variants pipeline as an sbatch job
# Processes *.fastq.gz files using the genome supplied.
# The genome file must end with the .fasta extension and was previously indexed via ./setup-escape-variants.sh.

ShowHelp()
{
   # Display Help
   echo "Runs a Slurm pipeline determining escape variants in fastq.gz files."
   echo
   echo "usage: $0 -g genome -i inputdir -m mode [-o outdir] [-w workdir] [-l logdir] [-e email] [-p project] [-j numjobs] [-d] [-D datetab]"
   echo "options:"
   echo "-g genome    *.fasta genome to use - required"
   echo "-i inputdir  directory containing *.fastq.gz files to process - required"
   echo "-m mode      run mode: 'h' (hospital surveillance) or 'c' (campus surveillance ) or 'e' (experimental) - required"
   echo "-o outdir    directory to hold output files - defaults to current directory"
   echo "-w workdir   directory that will hold a tempdir - defaults to current directory"
   echo "-l logdir    directory that will hold sbatch logs - defaults to /logs within outdir"
   echo "-e email     email address to notify on pipeline completion - defaults to empty(no email sent)"
   echo "-p project   name of the project (output filenames) - defaults to sars-cov2"
   echo "-j numjobs   number of array jobs to run in parallel - defaults to 4"
   echo "-D datetab   date.tab file used to create supermetadata.tab - defaults to skipping the supermetadata step"
   echo "-d           debug mode - skips deleting the tempdir"
   echo ""
   echo "NOTE: The input genome must first be indexed by running ./setup-variants-pipeline.sh."
   echo "NOTE: The inputdir, outdir, logdir, and workdir must be directories shared across the slurm cluster."
   echo ""
}

# set default argument values
export WORKDIR=$(pwd)
export OUTDIR=$(pwd)
export LOGDIR="$OUTDIR/logs"
export LOGSUFFIX=$$
export DELETE_EVTMPDIR=Y
export PROJECTNAME=sars-cov2
export EVMODE=""
export DATETAB=""
export SPIKEBED="spike.bed"

# parse arguments
while getopts "g:i:m:o:w:l:e:dp:j:D:" OPTION; do
    case $OPTION in
    g)
        export GENOME=$OPTARG
        ;;
    i)
        export INPUTDIR=$OPTARG
        ;;
    m)
        export EVMODE=$OPTARG
        ;;
    o)
        export OUTDIR=$OPTARG
        ;;
    w)
        export WORKDIR=$OPTARG
        ;;
    l)
        export LOGDIR=$OPTARG
        ;;
    e)
        export EMAIL=$OPTARG
        ;;
    p)
        export PROJECTNAME=$OPTARG
        ;;
    j)
        export MAX_ARRAY_JOBS=$OPTARG
        ;;
    d)
        export DELETE_EVTMPDIR=N
        ;;
    D)
        export DATETAB=$OPTARG
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
if [ -z "INPUTDIR" ]
then
   echo "ERROR: Missing required '-i inputdir' argument."
   echo ""
   ShowHelp
   exit 1
fi
# make sure mode is "h", "c" or "e"
if [[ "$EVMODE" != "h" && "$EVMODE" != "c" && "$EVMODE" != "e" ]]
then
   echo "ERROR: Required '-m mode' argument must be 'h', 'c' or 'e'."
   echo ""
   ShowHelp
   exit 1
fi

# check that the genome dictionary has been created by ./setup-escape-variants.sh
GENOME_BASE_NAME=$(basename $GENOME .fasta)
GENOME_DICTIONARY="$GENOME_BASE_NAME.dict"
if [[ ! -f "$GENOME_DICTIONARY" ]]
then
    echo "ERROR: Genome dictionary file $GENOME_DICTIONARY not found."
    echo "To fix run ./setup-escape-variants.sh -g $GENOME"
    echo ""
    exit 1
fi

# create logs directory to hold slurm logs (this directory must exist before sbatch can run)
mkdir -p $LOGDIR

# send email if user specifies to
SBATCH_FLAGS="$SBATCH_FLAGS --output=$LOGDIR/ev-pipeline-%j.out"
if [ ! -z "$EMAIL" ]
then
   echo "Emailing $EMAIL on pipeline completion."
   SBATCH_FLAGS="$SBATCH_FLAGS --mail-type=END --mail-user=$EMAIL"
fi

# run pipeline
JOBID=$(sbatch --parsable ${SBATCH_FLAGS} scripts/escape-variants-pipeline.sh)
echo "Submitted batch job $JOBID"
echo "To monitor main log run:"
echo "tail -f $LOGDIR/ev-pipeline-$JOBID.out"
