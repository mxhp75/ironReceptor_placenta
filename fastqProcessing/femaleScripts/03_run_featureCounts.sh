#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

# Notification configuration
#SBATCH --job-name=20231206_GDM_female_featureCounts
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20231206 ##

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

PROJROOT=/scratch/user/smit1924/20231206_GDM_female_grch38

##--------------------------------------------------------------------------------------------##
## Run samtools sort
##--------------------------------------------------------------------------------------------##

conda activate featureCounts # this is featureCounts v2.0.3

bash ${PROJROOT}/scripts/03_featureCounts.sh

conda deactivate
