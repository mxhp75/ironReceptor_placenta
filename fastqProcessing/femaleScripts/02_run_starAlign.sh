#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=16
#SBATCH --time=60:00:00
#SBATCH --mem=80GB

# Notification configuration
#SBATCH --job-name=20231206_GDM_female_starAlign
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20231206 ##

## Script for running STAR alignment ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load Miniconda3/4.9.2

##--------------------------------------------------------------------------------------------##
## Inside the 'STAR' conda environment
##--------------------------------------------------------------------------------------------##

#  + star   2.7.10b  h9ee0642_0  bioconda/linux-64
# openjdk-8.0.332            |       h166bdaf_0        97.8 MB  conda-forge
# picard-2.18.29             |                0        13.6 MB  bioconda

##--------------------------------------------------------------------------------------------##
## Assign variables
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/20231206_GDM_female_grch38

##--------------------------------------------------------------------------------------------##
## Run STAR alignment
##--------------------------------------------------------------------------------------------##

conda activate STAR

bash ${PROJROOT}/scripts/02_starAlign.sh

conda deactivate
