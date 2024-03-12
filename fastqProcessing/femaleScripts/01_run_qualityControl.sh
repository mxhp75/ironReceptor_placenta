#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=8
#SBATCH --time=40:00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --job-name=20240307_ironProjet_female_qualityControl
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240307 ##

## Script for running FastQC, fastp and MultiQC on raw fastq.gz RNA-seq data files ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load Miniconda3/4.9.2

##--------------------------------------------------------------------------------------------##
## Assign variables
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/ironReceptor_placenta/fastqProcessing

##--------------------------------------------------------------------------------------------##
## Run quality control script [FastQC (raw and trimmed files), fastp (adapter trim), MutiQC (combined report)
##--------------------------------------------------------------------------------------------##

conda activate qualityControl

bash ${PROJROOT}/femaleScripts/01_qualityControl.sh 

conda deactivate
