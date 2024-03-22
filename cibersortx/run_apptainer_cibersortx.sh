#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem=128GB

# Notification configuration
#SBATCH --job-name=smit1924_docker_cibersortx_ironReceptor
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240322 ##

## Script for running CIBERSORTx through apptainer -> docker ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load apptainer/1.1.9

##--------------------------------------------------------------------------------------------##
## Set environment variables
##--------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------##
## Set environment variables
##--------------------------------------------------------------------------------------------##

WORK_DIR="/scratch/user/smit1924/ironReceptor_placenta/cibersortx"
cd $WORK_DIR

##--------------------------------------------------------------------------------------------##
## Specify paths to data files and outputs
##--------------------------------------------------------------------------------------------##

OUT_DIR="${WORK_DIR}/outdir"

##--------------------------------------------------------------------------------------------##
## Run CIBERSORTx through apptainer
##--------------------------------------------------------------------------------------------##

apptainer run \
        -B ./data:/src/data \
        -B ${OUT_DIR}:/src/outdir \
        --pwd /src \
        fractions_latest.sif \
        --username melanie.smith@flinders.edu.au \
        --token e7d3e60a0c68bb01cbc06d543606066c \
        --verbose TRUE \
        --perm 100 \
        --fraction 0 \
        --single_cell TRUE \
        --refsample GSE182381_reference_sample.txt \
        --mixture tpm_bulkCountMatrix_all_gene_symbol.txt \
        --rmbatchSmode TRUE
