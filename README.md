# Functional consequences of reductive protein evolution in a minimal eukaryotic genome
## Authors: Jason Jiang, Rui Qu, Maria Grigorescu, Winnie Zhao, Aaron Reinke

Read the pre-print here: https://www.biorxiv.org/content/10.1101/2023.12.31.573788v2

This repo contains the environment and Snakemake pipeline needed to enact the main workflow in Figure 1.

## Prerequisites
- Local installation of Singularity >= 3.10
- Local installation of Python >= 3.10

## Set-up
### 1. Clone repo to local machine, then cd into repo

### 2. Run set-up script to initialize Python virtual environment for Snakemake

   ### NOTE: if you run into permission issues, try running chmod u+x setup.sh
   
   ./setup.sh
   
   source venv_snakemake/bin/activate

   snakemake --help

### 3. Initialize Singularity container for OrthoFinder, HMMER, cath-resolve-hits + all required R packages
   
   sudo singularity build container.sif container.def
   
### 4. Test Snakemake workflow + Singularity container
### NOTE: you may see a message like "System has not been booted with systemd as init system (PID 1). Can't operate."
### This does NOT affect our workflow, and only concerns datetime operations with R tidyverse (which we don't use)
   snakemake --cores all --use-singularity singularity_test
   

## Running the main workflow

## Ancilliary scripts to create figures (to be run interactively)
Note: these scripts require that you have already run the main workflow, so all orthogroups and domain architectures have already been assigned to your input proteomes.
