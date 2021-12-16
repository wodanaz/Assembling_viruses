# Escape Variants Pipeline

The escape variants pipeline is based on [Escape_Variants.md](https://github.com/wodanaz/Assembling_viruses/blob/main/Escape_Variants.md). 

## Requirements
- [bash](https://www.gnu.org/software/bash/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) version 6.4.1  

Snakemake 6.4.1 is required due to incompatible parameter changes in the latest version of snakemake.

## Installation
- Clone this repository
- Config File Setup
- Setup a Snakemake Environment
- Setup a DukeDSClient Environment - (Optional - used to stage data in and out)

### Config File Setup
A `config.sh` file must be created by copying the example file in place:
```
cp example-config.sh config.sh
```
This config file contains environment variables used by the pipeline.

`config.sh` Environment Variables
- `USE_MODULES` - Controls if the scripts use environment modules to provide software
  - Values "Y" and "N". Defaults to "Y"
- `FOREGROUND_MODE` - Runs snakemake in the foreground.
  - Values "Y" and "N".  Defaults to "N"
  - This feature is for Slurm clusters that do not support `sacct` on worker nodes
- `SM_CONDA_PREFIX` - Environment variable to overide the conda cache directory.
  - Defaults to `conda` within the data directory.

### Snakemake Environment
Snakemake must be installed in a conda environment named `snakemake`.

To create this environment using the `conda` command run:
```
conda create -n snakemake -c conda-forge -c bioconda mamba snakemake==6.4.1
```

On the HARDAC cluster the the `snakemake` environment can be created on an interactive node like so:
```
module load Anaconda3/2019.10-gcb02
conda create -n snakemake -c conda-forge -c bioconda mamba snakemake==6.4.1
```

### DukeDSClient Environment
If not using environment modules a conda environment named `dds` should be created by running:
```
conda env create -f dds-environment.yml
```
When using environment modules (`USE_MODULES=Y`) the `ddsclient` environment module will be used.


## Running the pipeline
The pipeline uses a `data` directory to organize input, conda environments, and output.
By default the pipeline stages data in from DDS so a project name is required.
The main pipeline script has two required arguments:
```
./run-escape-variants.sh -d <datadir> -i <projectname>
```
So to use "data" as the name of your `data` directory and "CVTestSamples" as your project name run:
```
./run-escape-variants.sh -d data -i CVTestSamples
```
The above command will:
- Download data into `data/input/CVTestSamples`
- Run snakemake storing output files in `data/output/CVTestSamples_snakemake`
- Copy output files from `data/output/CVTestSamples_snakemake` to `data/output/CVTestSamples`
- Upload results into a project named `CVTestSamples_results`

If a problem occurs when running the pipeline after making a change (such as increasing memory requested) the same `run-escape-variants.sh` command can be run again.
Snakemake will restart where it left off.

### Skipping staging data
If you do not want to download and upload data to DDS put the input fastq.gz files in the
appropriate `<datadir>/input/<projectname>` directory and add the `-s -S` arguments to `run-escape-variants.sh`.
- `-s` disables staging data in (downloading) from DDS
- `-S` disables staging data out (upload) to DDS

For example using data as `datadir` and `CVTestSamples` as a project name input files should be placed in `data/input/CVTestSamples` and the following command run:
```
./run-escape-variants.sh -d data -i CVTestSamples -s -S
```


### Email on pipeline finishing
The `-e <email>` argument will cause the `run-escape-variants.sh` script to send an email when the job completes.

For example:
```
./run-escape-variants.sh -d data -i CVTestSamples -e yourusername@example.com
```
  
### Preserve Input and Intermediate Data
By default the pipeline will delete the input directory and intermediate files/directories. 
If you want to preserve the input and intermediate directories pass the `-k` argument.

### Run Script Help
You can see command line help by running `./run-escape-variants.sh` without arguments.

## Running Snakemake Directly
The snakemake pipeline can be run directly bypassing the helper scripts used above.
This can be helpful when making changes to the pipeline.
  
### Configure the Pipeline
All settings for the pipeline are stored in the [config/config.yaml](config/config.yaml) file.

Configuration options:
- __project__ - Provides a unique name to be used when creating output files
- __genome__ - File path to a genome FASTA file used in the pipeline. Default value `resources/NC_045512.fasta`.
- __spike__ -  File path to a bed file containing the spike region. Default value `resources/spike.bed`.
- __inputdir__ - directory containing input fastq.gz sequences to process. Default value `""` (the base directory of this repo).
If __inputdir__ is changed to a path other than `""` the genome, spike, datetab should be set to absolute paths.
Otherwise the workflow will not be able to find these files.

### Run Snakemake
Copy your *.fastq.gz samples into the base directory of this repo.
From within the base directory of this repo run snakemake:
```
snakemake --cores 1 --use-conda
```
The above command will use one core and create conda environments for the tools used in the pipeline.
The number of cores can be increased to improve throughput.
After the command finishes the pipeline output files and directories will be in the base directory of this repo.

#### Running snakemake with Slurm
The snakemake command creates multiple jobs and monitors their progress.
Since the snakemake command runs in the foreground it is best to create an sbatch script to run it.

For example to run with 10 concurrent processes, using conda for the required tools, and the slurm profile you
could create a sbatch script named `ev.sh` with the following content:
```
#!/bin/bash
module load Anaconda3/2019.10-gcb02
conda activate snakemake
snakemake --cores 10 --use-conda --profile slurm/
```
Then run the sbatch script like so:
```
sbatch ev.sh
```
If the environment module is not available on your cluster remove the `module load Anaconda3/2019.10-gcb02` line and setup conda.


## Using local storage on a single node
The pipeline can be run on single Slurm node utilizing local storage.
Using this mode may help if a Slurm cluster's shared storage is a bottleneck.
When using this mode snakemake, conda environments, input data, and intermediate data will all be stored on local storage.
Only the output files and top level logs will stored on the shared storage.
The top level log files must be on shared storage for the inital sbatch job to run.

### Setup config
Running in the single node mode requires special configuration.
To create the config files run `./setup-config.sh` and answer the prompts.
This script creates/replaces the contents of the `config.sh` and `smprofile/config.yaml` files.
The script will prompt for the following items:
- nodename - the node where all jobs will be run
- scratch base directory - directory where , account, partition
- account - Slurm account to use for the jobs
- partition - Slurm partition to use for the jobs (must include the node specified above)


### Setup miniconda3
Miniconda3 must be installed on the local(scratch) storage on the node to run the workflow.
The Snakemake and DukeDSClient environments installed in the appropriate environments within the miniconda3 installation.

### Running on a single node
To run on a single node requires specifing the optional `-o` output directory argument.
This argument specifies the location for output files.
So if `/scratch/user1/data` is the local storage path and `/hpc/group/mygroup/user1` is the shared output files path.
The pipeline could be run on project `CVTestSamples` with 120 jobs like so:
```
./run-escape-variants.sh -d /scratch/user1/data -o /hpc/group/mygroup/user1 -i CVTestSamples -j 120
```

