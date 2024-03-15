#!/bin/bash -l

## Script to run STAR with gene and transcript level counts ##
## Requires STAR index to be pre-prepared ##
## Melanie Smith ##
## 20240313 ##

##--------------------------------------------------------------------------------------------##
## set project root directory
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/ironReceptor_placenta/fastqProcessing

##--------------------------------------------------------------------------------------------##
## set project sub-directories
##--------------------------------------------------------------------------------------------##

mkdir -p ${PROJROOT}/aligned_data

refs=/scratch/user/smit1924/refSeq # all the reference datasets are here
TRIMFQ=${PROJROOT}/femaleOutput/trim_data    # the trimmed fastq files are here
ALIGN_OUTPUT=${PROJROOT}/femaleOutput/aligned_data	# put the aligned data here

##--------------------------------------------------------------------------------------------##
## create align_data sub-directory
##--------------------------------------------------------------------------------------------##

mkdir -p ${PROJROOT}/femaleOutput/aligned_data

##--------------------------------------------------------------------------------------------##
## Check the aligned data directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${ALIGN_OUTPUT} ]]
  then
    echo -e "${ALIGN_OUTPUT} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${ALIGN_OUTPUT} \n"

##--------------------------------------------------------------------------------------------##
## Set a few compute specifications
##--------------------------------------------------------------------------------------------##

PICARD=$EBROOTPICARD/picard.jar

build=GRCh38
cores=16

##--------------------------------------------------------------------------------------------##
## Align trimmed data to reference genome
##--------------------------------------------------------------------------------------------##

for FQGZ in ${TRIMFQ}/*_R1_001_trim.fastq.gz; do

 SampleName=$(basename ${FQGZ} _R1_001_trim.fastq.gz)

        STAR --genomeDir ${refs}/STAR_index/GRCh38_female \
             --readFilesIn ${FQGZ} ${FQGZ/R1/R2} \
             --readFilesCommand zcat \
             --quantMode GeneCounts \
	     --genomeLoad NoSharedMemory \
	     --outFilterType BySJout \
             --outFilterMismatchNmax 999 \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix ${ALIGN_OUTPUT}/"${SampleName}"_"${build}"_ \
             --outSAMattrRGline ID:"${SampleName}" \
                                LB:library \
                                PL:illumina \
                                PU:machine \
                                SM:"${build}" \
             --outSAMmapqUnique 60 \
             --runThreadN ${cores}
done
