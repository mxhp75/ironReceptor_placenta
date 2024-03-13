#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=8
#SBATCH --time=40:00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --job-name=20240307_ironProjet_female_truncation_test
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240313 ##

## Script for running a quick truncation test on the raw fastq.gz RNA-seq data files ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared

##--------------------------------------------------------------------------------------------##
## Assign variables
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/ironReceptor_placenta/fastqProcessing
output_file=${PROJROOT}/femaleOutput/rawQC/truncation_test.txt
output_file_error=${PROJROOT}/femaleOutput/rawQC/truncation_errorOnly.txt
input_directory=${PROJROOT}/female_fastq

##--------------------------------------------------------------------------------------------##
## Run truncation test on all fastq.gz files in the input directory and output short report
##--------------------------------------------------------------------------------------------##

# Iterate over files
#for file in ${input_directory}/SAGCFN_22_*.fastq.gz; do
#    echo "Checking file: $file" >> "$output_file"
#    if zcat "$file" | tail >> "$output_file" 2>&1; then
#        echo "-------------------------------------" >> "$output_file"
#    else
#        echo "-------------------------------------" >> "$output_file"
#        echo "Error: Unexpected end of file in $file" >> "$output_file"
#    fi
#done

# Iterate over files -> error only
for file in ${input_directory}/SAGCFN_22_*.fastq.gz; do
    echo "Checking file: $file" >> "$output_file_error"
    if ! zcat "$file" >/dev/null 2>&1; then
        echo "$file: Error - Unexpected end of file" >> "$output_file_error"
    fi
done


