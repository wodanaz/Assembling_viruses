# Escape Variants Pipeline

The escape variants pipeline is based on [Escape_Variants.md](https://github.com/wodanaz/Assembling_viruses/blob/main/Escape_Variants.md). 

## Requirements
- [bash](https://www.gnu.org/software/bash/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)  

## Installation
Clone this repository.

## Configure the Pipeline
All settings for the pipeline are stored in the [config/config.yaml](config/config.yaml) file.

Configuration options:
- __project__ - Provides a unique name to be used when creating output files
- __genome__ - File path to a genome FASTA file used in the pipeline. Default value `resources/NC_045512.fasta`.
- __spike__ -  File path to a bed file containing the spike region. Default value `resources/spike.bed`.
- __inputdir__ - directory containing input fastq.gz sequences to process. Default value `""` (the base directory of this repo).
- __mode__ - Determines which version of the pipeline to run. Must be "campus", "hospital", or "experimental". Default value "campus".
- __datetab__ - Optional path to a tsv file with two columns: SampleName and Date. Not default value.
If __inputdir__ is changed to a path other than `""` the genome, spike, datetab should be set to absolute paths.
Otherwise the workflow will not be able to find these files.

## Run the Pipeline
Copy your *.fastq.gz samples into the base directory of this repo.
From within the base directory of this repo run snakemake:
```
snakemake --cores 1 --use-conda
```
The above command will use one core and create conda environments for the tools used in the pipeline.
The number of cores can be increased to improve throughput.
After the command finishes the pipeline output files and directories will be in the base directory of this repo.

### HARDAC/Slurm Instructions
One way to scale the pipeline beyond a single machine is to use a Slurm cluster.

#### Create a conda environment for snakemake
On the HARDAC cluster the Anaconda3 module provides the conda command.
The conda command can be used to create an environment containing the snakemake tool.
From a HARDAC session run the following script to create a __snakemake__ conda environment:
```
srun scripts/create-snakemake-env.sh
```

#### Running snakemake with sbatch
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

#### Staging data with Duke Data Service
The `run-dds-escape-variants.sh` script stages data to/from projects in the [Duke Data Service](https://dataservice.duke.edu/). This script will download a DDS project, run the escape variants pipeline on any fastq.gz files, then uploads the results to a new DDS project. The resulting project name is created by adding "_results" to the input project name.

##### Additional requirements
The `ddsclient/3.2.0-gcb01` environment module is required for uploading and downloading from Duke Data Service.
DDS credentials need to be placed in a `~/.ddsclient` file.
See the [documentation for creating the DukeDSClient config file](https://github.com/Duke-GCB/DukeDSClient#config-file-setup).

##### Data Directory
The `run-dds-escape-variants.sh` script requires a `datadir` argument. `run-dds-escape-variants.sh` will create two subdirectories in `datadir` named `input` and `output`. The `input` directory will contain a subdirectory for each project downloaded from DDS. These project specific input directories will be deleted when the pipeline successfully completes. The `output` directory will also contain a project specific directory containing the output files from a pipeline and the associated logs. The `output` project specific subdirectory is not deleted.

##### Example
Given the following for your current directory:
1. A directory used to hold input and output data named `datadir`
2. An indexed genome `MT246667.fasta` file
3. The name of a DDS project to download in this case `sars-cov2-example`.
```
sbatch run-dds-escape-variants.sh -d datadir -g MT246667.fasta -i sars-cov2-example
```
The surveillance mode can be activated by passing the `-s` flag.

##### Email on Completion
To receive emails when the top level sbatch job completes pass the standard sbatch email flags.

Example:
```
sbatch --mail-type=END --mail-user=<email> run-dds-escape-variants.sh -d datadir -g MT246667.fasta -i sars-cov2-example
```

##### Help
You can see command line help by running `./run-dds-escape-variants.sh` without arguments.
