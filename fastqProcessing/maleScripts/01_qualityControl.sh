#!/bin/bash -l

## Melanie Smith ##
## 20240307 ##

## Script for performing FastQC, fastp and MultiQC on raw fastq.gz RNA-seq data files ##
## To be run on DeepThought ##

# nb: conda environment activated in the run script

##--------------------------------------------------------------------------------------------##
## set project root directory
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/ironReceptor_placenta/fastqProcessing

##--------------------------------------------------------------------------------------------##
## set up log file
##--------------------------------------------------------------------------------------------##

mkdir -p ${PROJROOT}/maleOutput/log
exec 1> ${PROJROOT}/maleOutput/log/qualitycontrol.log 2>&1
set -ex

##--------------------------------------------------------------------------------------------##
## create qc and trim data sub-directories
##--------------------------------------------------------------------------------------------##

mkdir -p ${PROJROOT}/maleOutput/trim_data/ ${PROJROOT}/maleOutput/trimQC/multiqc ${PROJROOT}/maleOutput/rawQC/multiqc

##--------------------------------------------------------------------------------------------##
## set number of threads
##--------------------------------------------------------------------------------------------##

THREADS=8

##--------------------------------------------------------------------------------------------##
## set variables for project sub-directories
##--------------------------------------------------------------------------------------------##

RAWFQ=${PROJROOT}/male_fastq		# this is where the untrimmed, merged fastq files are
RAWQC=${PROJROOT}/maleOutput/rawQC	# put the raw QC reports here
TRIMFQ=${PROJROOT}/maleOutput/trim_data	# put the trimmed fastq files here
TRIMQC=${PROJROOT}/maleOutput/trimQC	# put the trimmed QC reports here

##--------------------------------------------------------------------------------------------##
## Check the project root directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${PROJROOT} ]]
  then
    echo -e "${PROJROOT} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${PROJROOT}\n"


##--------------------------------------------------------------------------------------------##
## Check the raw directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${RAWFQ} ]]
  then
    echo -e "${RAWFQ} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${RAWFQ}\n"


##--------------------------------------------------------------------------------------------##
## Check the raw data QC directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${RAWQC} ]]
  then
    echo -e "${RAWQC} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${RAWQC}\n"


##--------------------------------------------------------------------------------------------##
## Check the trimmed data directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${TRIMFQ} ]]
  then
    echo -e "${TRIMFQ} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${TRIMFQ}\n"


##--------------------------------------------------------------------------------------------##
## Check the trimmed data QC directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${TRIMQC} ]]
  then
    echo -e "${TRIMQC} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${TRIMQC}\n"


##--------------------------------------------------------------------------------------------##
## run FastQC on raw data, unless already preformed (change file names as needed)
##--------------------------------------------------------------------------------------------##

for R1 in ${RAWFQ}/*_R1_001.fastq.gz; do
  echo -e "Found ${R1}\n"
  R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz
  echo -e "The R2 file should be ${R2}\n"
  base1=$(basename $R1)
  base2=$(basename $R2)

  for F in $R1 $R2; do
    echo $F
    if [[ ! -f ${RAWQC}/${baseF%.fastq.gz}_fastqc.html  ]]
      then
        echo -e "Running FastQC on ${F}\n"
        fastqc -o ${RAWQC} -t ${THREADS} ${F}
      else
        echo -e "FastQC already performed on ${F}, skipping\n"
    fi
  done


##--------------------------------------------------------------------------------------------##
## trim reads to remove adapter sequence and discard poor quality reads
##--------------------------------------------------------------------------------------------##

  if [[ ! -e ${TRIMFQ}/${base1%.fastq.gz}_trim.fastq.gz || ! -e ${TRIMFQ}/${base2%.fastq.gz}_trim.fastq.gz ]]
    then
      echo -e "Running fastp on ${F}\n"
      fastp \
        -i ${R1} -I ${R2} \
        -o ${TRIMFQ}/${base1%.fastq.gz}_trim.fastq.gz \
        -O ${TRIMFQ}/${base2%.fastq.gz}_trim.fastq.gz \
        --cut_right --cut_window_size 4 --cut_mean_quality 20 \
        --length_required 75
  fi
done


##--------------------------------------------------------------------------------------------##
## find forward and reverse trimmed read files (change file names as needed)
##--------------------------------------------------------------------------------------------##

for TR1 in ${TRIMFQ}/*_R1_001_trim.fastq.gz; do
  echo -e "Found ${TR1}\n"
  TR2=${TR1%_R1_001_trim.fastq.gz}_R2_001_trim.fastq.gz
  echo -e "The TR2 file should be ${TR2}\n"

##--------------------------------------------------------------------------------------------##
## run FastQC on trimmed reads, unless already performed
##--------------------------------------------------------------------------------------------##

  for F in $TR1 $TR2; do
    echo $F
    baseF=$(basename $F)
    if [[ ! -f ${TRIMQC}/${baseF%.fastq.gz}_fastqc.html  ]]
      then
        echo -e "Running FastQC on ${baseF}\n"
        fastqc -o ${TRIMQC} -t ${THREADS} ${F}
      else
        echo -e "FastQC already performed on ${F}, skipping\n"
    fi
  done
done

##--------------------------------------------------------------------------------------------##
## collate raw QC data into multiQC report
##--------------------------------------------------------------------------------------------##

if [[ ! -f ${RAWQC}/multiqc/multiqc-report_rawdata.html  ]]
  then
    multiqc \
    --outdir ${RAWQC}/multiqc \
    --filename multiqc-report_rawdata.html \
    ${RAWQC}
fi

##--------------------------------------------------------------------------------------------##
## collate trimmed data into multiQC report
##--------------------------------------------------------------------------------------------##

if [[ ! -f ${TRIMQC}/multiqc/multiqc-report_trimdata.html  ]]
  then
    multiqc \
    --outdir ${TRIMQC}/multiqc \
    --filename multiqc-report_trimdata.html \
    ${TRIMQC}
fi
