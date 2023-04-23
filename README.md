# map_fastqs_snakemake

A program for mapping fastq files to a reference genome with bwa-mem, then
converting the output sam file to a sorted indexed bam file.

This output file is suitable for viewing in UCSC genome browser or IGV, or for
passing on to other analysis pipelines that need sorted indexed bam files.

This pipeline uses conda environments to standardize the software used for each
run.

## Installation:
  - Install conda: https://github.com/conda-forge/miniforge#mambaforge
  - Create a conda environment and install snakemake there:
  ```bash
  conda create -c conda-forge -c bioconda -n snakemake snakemake
  conda activate snakemake
  ```
  - setup the conda environment to use strict mode:
  ```bash
  conda config --set channel_priority strict
  ```

### Setup your environment:
  - Change directory to a folder where you want to run the analysis
  - clone this git repository into the folder

## Usage:
  - Edit the config.yaml file using the instructions in the comments. Use a text
  editor that outputs unix line endings (e.g. vscode, notepad++, gedit, micro,
  emacs, vim, vi, etc.)
  - If snakemake is not your active conda environment, activate snakemake with:
  ```bash
  conda activate snakemake
  ```
  - Run snakemake (if you have mamba) with:
  ```bash
  snakemake -s map_fastqs.smk --cores 4 --use-conda
  ```
  - or (if you're using some other conda) with:
  ```bash
  snakemake -s map_fastqs.smk --cores 4 --use-conda --conda-frontend conda
  ```
