!/bin/bash

# Oct 2022
# Carmel Maher
# This script is to submit a FastQC analysis on trimmed reads through Flinders' DeepThought HPC Slurm
# This script therefore assumes submission to Flinders Deepthought with an assumed directory structure and use of available modules


#------------------------------------
# Required Modules:

#Required Installations
#fastQC (v0.11.9)
#openjdk::conda-forge/linux-64::openjdk-17.0.3-h85293d2_2


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

cd /scratch/user/mahe0050/DE-analysis/1_trimmedData/q/ || error_exit "$LINENO: directory error 1"

for file in *R1_cleanq.fq.gz
do
	echo $file
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo "FILESTEM2" $FILESTEM2
	
	fastqc --noextract --threads 8 -o ./FastQC_Cleanq ./$FILESTEM2"_R1_cleanq.fq.gz" ./$FILESTEM2"_R2_cleanq.fq.gz" || error_exit "$LINENO: Error with FastQC-q at "$FILESTEM2""
done

	
cd /scratch/user/mahe0050/DE-analysis/1_trimmedData/ || error_exit "$LINENO: directory error 2"

for file in *R1_clean.fq.gz
do
	echo $file
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo "FILESTEM2" $FILESTEM2
	
	fastqc --noextract --threads 8 -o ./FastQC_Clean ./$FILESTEM2"_R1_clean.fq.gz" ./$FILESTEM2"_R2_clean.fq.gz" || error_exit "$LINENO: Error with FastQC at "$FILESTEM2""
done