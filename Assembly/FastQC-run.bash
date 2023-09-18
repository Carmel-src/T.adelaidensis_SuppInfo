!/bin/bash

# November 2019
# Carmel Maher
# This script is to submit a FastQC analysis on raw reads through Flinders' DeepThought HPC Slurm
# This script therefore assumes submission to Flinders Deepthought with an assumed directory structure and use of available modules


#------------------------------------
# Required Modules:
# module add fastqc/0.11.8


#------------------------------------
# Terry's error exit

function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

#------------------------------------

cd /scratch/user/mahe0050/DE-analysis/0_rawData/AGRF_CAGRF13871_HCGYYBCXY || error_exit "$LINENO: directory error 1"
for file in *.fastq.gz
do
	FILESTEM=${file%.*}
	echo $FILESTEM

	fastqc $FILESTEM.gz --extract --outdir /scratch/user/mahe0050/DE-analysis/0_rawData/FastQC || error_exit "$LINENO: Error with FastQC at "$FILESTEM""
done
	
cd /scratch/user/mahe0050/DE-analysis/0_rawData/AGRF_CAGRF20022_CDNJBANXX || error_exit "$LINENO: directory error 2"
for file in *.fastq.gz
do
	FILESTEM=${file%.*}
	echo $FILESTEM

	fastqc $FILESTEM.gz --extract --outdir /scratch/user/mahe0050/DE-analysis/0_rawData/FastQC || error_exit "$LINENO: Error with FastQC at "$FILESTEM""
done
