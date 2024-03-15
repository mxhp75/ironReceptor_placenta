#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

# Notification configuration
#SBATCH --job-name=20240315_ironProject_male_featureCounts
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240315 ##

## Script for running Feature Counts ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load Miniconda3/4.9.2

##--------------------------------------------------------------------------------------------##
## Inside the 'featureCounts' conda environment
##--------------------------------------------------------------------------------------------##

# mamba
# subread
# multiqc
# samtools

##--------------------------------------------------------------------------------------------##
## Assigni project root variables
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/ironReceptor_placenta

##--------------------------------------------------------------------------------------------##
## Run samtools sort
##--------------------------------------------------------------------------------------------##

conda activate featureCounts # this is featureCounts v2.0.3

bash ${PROJROOT}/fastqProcessing/maleScripts/03_featureCounts.sh

conda deactivate
