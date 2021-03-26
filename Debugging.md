
# Debugging
The [scripts/escape-variants-pipeline.sh](https://github.com/wodanaz/Assembling_viruses/blob/main/scripts/escape-variants-pipeline.sh) script manages the escape variants pipeline by running a series of sbatch scripts. If one of these sbatch script fails the script can be run directly. An interactive session and some environment variables must be setup before directly running these scripts.

## Create an Interactive Session
Start an interactive job with enough memory for the script you want to run.

For example:
```
srun -p interactive --pty --mem 10G bash
```

## Activate the Conda Environment
```
module load Anaconda3/2019.10-gcb02
conda activate escapevariants
```

## General Environment Variables
The scripts require environment variables to be exported before running them.

The following environment variables may need to set:
- __EVDIR__ - directory used to hold intermediate and output files
- __GENOME__ - path to indexed FASTA genome file for scripts that use it

For example to use `tmp.PnkjJwoLPx` as an intermediate/output directory and `NC_045512.fasta` as the genome file run the following:
```
export EVDIR=tmp.PnkjJwoLPx
export GENOME=NC_045512.fasta
```

## Script Types
There are two types of sbatch scripts called from [scripts/escape-variants-pipeline.sh](https://github.com/wodanaz/Assembling_viruses/blob/main/scripts/escape-variants-pipeline.sh):
- __simple scripts__ that do not have any command line arguments
- __array job scripts__ that have a single list file command line argument

You can determine the type of script by seeing how it is called from scripts/escape-variants-pipeline.sh.

The __simple scripts__ are called using sbatch directly like so:
```
sbatch --wait "--output=${LOGDIR}/create-coverage-table-%j.out" \
    ./scripts/create-coverage-table.sh
```
The __array job scripts__ are called using the sbatch-array.sh helper script:
```
./scripts/sbatch-array.sh \
    ./scripts/remove-nextera-adapters.sh $EVDIR/reads.list
```

## Run Simple Script
The __simple scripts__ can be run directly without any arguments.
For example `create-coverage-table.sh` is called in `escape-variants-pipeline.sh` like so:
```
sbatch --wait "--output=${LOGDIR}/create-coverage-table-%j.out" \
    ./scripts/create-coverage-table.sh
```
To manually run this script without sbatch run the following command:
```
./scripts/create-coverage-table.sh
```

## Run Array Job Script
__Array job scripts__ receive a single command line argument of a ".list" file.
These ".list" files contain a list of filenames to be processed by the __array job script__.
Normally the __array job script__ is run as an sbatch array job processing multiple files in parallel.
To run these scripts directly an environment variable must be setup specifiying which filename in the ".list" file to process.

### Specify which file to process
The __SLURM_ARRAY_TASK_ID__ environment variable controls which file is processed by a __array job script__.
Set the __SLURM_ARRAY_TASK_ID__ environment variable to an integer specifying which file in the list to process.
For example to enable processing for the first filename in the ".list" file run the following:
```
export SLURM_ARRAY_TASK_ID=1
```

### Run the Script
The argument to pass to a __array job script__ can be determined by looking at how the script is called in `escape-variants-pipeline.sh`.

For example `remove-nextera-adapters.sh` is called from `escape-variants-pipeline.sh` like so:
```
./scripts/sbatch-array.sh \
    ./scripts/remove-nextera-adapters.sh $EVDIR/reads.list
```
To manually run this script remove the `sbatch-array.sh` part and run it like so:
```
./scripts/remove-nextera-adapters.sh $EVDIR/reads.list
```

## Use ev-pipeline log as a guide
When the escape variants pipeline is run normally a log file named ev-pipeline*.log will be created.
This log includes indented commands that were run. 
After creating an interactive session and setting up __SLURM_ARRAY_TASK_ID__ you can paste the indented commands from this log into your terminal.

Here is an excerpt from one of these logs:
```
...

Current directory is /data/itlab/Assembling_viruses

Creating a temp directory
Running using directory /data/itlab/Assembling_viruses/tmp.PnkjJwoLPx

Exported Environment Variables
    export EVDIR=/data/itlab/Assembling_viruses/tmp.PnkjJwoLPx
    export GENOME=NC_045512.fasta
    export INPUTDIR=.

...

GATK Step 6b How much of the reference genome is covered by more than 1 read?
Processing 1 file(s) using:
    ./scripts/run-depth-coverage.sh /data/itlab/Assembling_viruses/tmp.PnkjJwoLPx/bqsrs.list
Submitted batch job 24786344
GATK Step 6b - Done

...
```


1. Create an interactive session and set up __SLURM_ARRAY_TASK_ID__
```
srun -p interactive --pty --mem 10G bash
export SLURM_ARRAY_TASK_ID=1
```

2. Change to the directory of this repo if necessary
```
cd /data/itlab/Assembling_viruses
```

3. Export the environment variables by running the commands from the log
```
export EVDIR=/data/itlab/Assembling_viruses/tmp.PnkjJwoLPx
export GENOME=NC_045512.fasta
export INPUTDIR=.
```

4. Run the command with the appropriate input file based on the log
```
./scripts/run-depth-coverage.sh /data/itlab/Assembling_viruses/tmp.PnkjJwoLPx/bqsrs.list
```
