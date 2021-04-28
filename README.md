# Escape Variants Slurm Pipeline

The escape variants Slurm pipeline is based on [Escape_Variants.md](https://github.com/wodanaz/Assembling_viruses/blob/main/Escape_Variants.md). 

## Requirements
- [bash](https://www.gnu.org/software/bash/)
- [Slurm cluster](https://slurm.schedmd.com/) with shared storage
- [conda](https://conda.io/projects/conda/en/latest/index.html)  

_By default __conda__ is provided via the [Anaconda3/2019.10-gcb02](https://github.com/Duke-GCB/helmod/blob/master/rpmbuild/SPECS/Anaconda3-2019.10-gcb02.spec) environment module. This module name can be overridden via the `ANACONDAMODULE` environment variable._

### Conda Environment
The pipeline requires a conda environment named `escapevariants` containing the following software:
  - seqtk=1.3
  - pangolin=2.3.2
  - gatk4=4.1.3.0
  - trim-galore=0.6.5
  - bcftools=1.10.2
  - bedtools=2.25.0
  - bwa=0.7.12
  - perl=5.32.0
  - cutadapt=2.3
  - picard=2.25.0
  - samtools=1.10
  - tabix=0.2.6

## Installation
Clone this repository onto a shared location in your Slurm cluster. 
Create the config.sh file containing global configuration settings such as the slurm account.
Then create the `escapevariants` conda environment.

### Create config.sh
The config.sh file can be created by copying the example file in place.
```
cp example-config.sh config.sh
```

### HARDAC Installation Instructions
On HARDAC the Anaconda3 module provides the conda command.
From an interactive session you can use the Anaconda3 module to create the `escapevariants` conda environment by running:
```
module load Anaconda3/2019.10-gcb02
conda env create -f environment.yml
```

### General Installation Instructions
Using an Anaconda or Miniconda installallation create the `escapevariants` conda environment by running:
```
conda env create -f environment.yml
```

## Setup
Before the pipeline can be run the input genome FASTA file must be processed by [setup-escape-variants.sh](https://github.com/wodanaz/Assembling_viruses/blob/main/setup-escape-variants.sh).
The setup script will index the reference genome and create a dictionary file used by the escape variants pipeline. The output files will be created in the same directory as the FASTA genome file.
For example to process a genome file named `MT246667.fasta` run:
```
./setup-escape-variants.sh -g MT246667.fasta
```

## Running
The `run-escape-variants.sh` script is used to run the pipeline on fastq.gz files already staged onto your cluster.
The script should be passed arguments specifying:
1. a directory containing input fastq.gz files
2. an indexed genome file
3. an output directory

The script submits a main sbatch job that will run additional sbatch jobs for each part of the pipeline. A temporary directory will be created to hold intermediate files. The temporary directory will be deleted when the pipeline successfully completes.

### Example
Given the following:
1. input fastq.gz files are stored in a directory named `inputdir`
2. an indexed genome `MT246667.fasta` file
3. output files will be written to a directory named `results` 

Run the pipeline like so:
```
./run-escape-variants.sh -g MT246667.fasta -i inputdir -o results
```

#### Surveillance Mode
By default the pipeline compiles all tab tables into one for depth and genotype.
This step can be skipped by turning on surveillance mode by passing the `-s` flag.

Example:
```
./run-escape-variants.sh -g MT246667.fasta -i inputdir -o results -s
```

#### Email on Completion
To receive emails when the main sbatch job completes pass your email after the `-e` flag.

Example:
```
./run-escape-variants.sh -g MT246667.fasta -i inputdir -o results -e <email>
```

#### Debugging
The pipeline creates a temp directory that contains all intermediate files. By default this directory is deleted when the pipeline completes successfully.
If you want to preserve the temp directory pass the `-d` debug argument.

Example:
```
./run-escape-variants.sh -g MT246667.fasta -i inputdir -o results -d
```

#### Help
To see command line help by running `./run-escape-variants.sh` without arguments.
There are additional command line arguments for controlling where the logs and temp directory are created.


## Staging data with Duke Data Service 
The `run-dds-escape-variants.sh` script stages data to/from projects in the [Duke Data Service](https://dataservice.duke.edu/). This script will download a DDS project, run the escape variants pipeline on any fastq.gz files, then uploads the results to a new DDS project. The resulting project name is created by adding "_results" to the input project name.

### Additional requirements
The `ddsclient/3.2.0-gcb01` environment module is required for uploading and downloading from Duke Data Service.
DDS credentials need to be placed in a `~/.ddsclient` file.
See the [documentation for creating the DukeDSClient config file](https://github.com/Duke-GCB/DukeDSClient#config-file-setup).

### Data Directory
The `run-dds-escape-variants.sh` script requires a `datadir` argument. `run-dds-escape-variants.sh` will create two subdirectories in `datadir` named `input` and `output`. The `input` directory will contain a subdirectory for each project downloaded from DDS. These project specific input directories will be deleted when the pipeline successfully completes. The `output` directory will also contain a project specific directory containing the output files from a pipeline and the associated logs. The `output` project specific subdirectory is not deleted.

### Example
Given the following for your current directory:
1. A directory used to hold input and output data named `datadir`
2. An indexed genome `MT246667.fasta` file
3. The name of a DDS project to download in this case `sars-cov2-example`.
```
sbatch run-dds-escape-variants.sh -d datadir -g MT246667.fasta -i sars-cov2-example
```
The surveillance mode can be activated by passing the `-s` flag.

#### Email on Completion
To receive emails when the top level sbatch job completes pass the standard sbatch email flags.

Example:
```
sbatch --mail-type=END --mail-user=<email> run-dds-escape-variants.sh -d datadir -g MT246667.fasta -i sars-cov2-example
```


#### Help
You can see command line help by running `./run-dds-escape-variants.sh` without arguments.

