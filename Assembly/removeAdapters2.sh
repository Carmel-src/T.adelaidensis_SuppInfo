#!/bin/bash

# October 2022
# Carmel Maher
# This script is to submit a Cutadapt script to cut adapters and quality on raw reads through Flinders' DeepThought HPC Slurm
# This script therefore assumes submission to Flinders Deepthought with an assumed directory structure and use of available modules


#------------------------------------
#Run from within ShortConda environment 
#module load Miniconda3/4.9.2
#conda activate ShortConda (this is input in every SLURM script)

#Required Installations
#cutadapt (version 4.1)
#fastQC (v0.11.9)

#------------------------------------

# Illumina Truseq Primers:
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA (to remove from R1)
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT (to remove from R2)

# Clontech SMARTer Primers:
# SMARTer II A Oligonucleotide
# 5'– AAGCAGTGGTATCAACGCAGAGTACXXXXX –3' (X = undisclosed base in the proprietary SMARTer oligo sequence)
# 3' SMART CDS Primer II A
# 5’– AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N –3’ (N = A, C, G, or T; N-1 = A, G, or C)
# 
# FWDPRIMER 	= 	AAGCAGTGGTATCAACGCAGAGTACNNNNN
# RCFWDPRIMER 	= 	NNNNNGTACTCTGCGTTGATACCACTGCTT
# REVPRIMER 	= 	AAGCAGTGGTATCAACGCAGAGTACT
# RCREVPRIMER 	= 	AGTACTCTGCGTTGATACCACTGCTT


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

#	in = /scratch/user/mahe0050/DE-analysis/0_rawData
#	out = /scratch/user/mahe0050/DE-analysis/1_trimmedData

#------------------------------------

#this script assumes you start IN the clean data folder as per DeepThought's paralell file system
#cd $BGFS


cd /scratch/user/mahe0050/DE-analysis/0_rawData || error_exit "$LINENO: directory error 1"


for file in *R1.fastq.gz
do
	echo $file
	FILESTEM=${file%_*}
	echo $FILESTEM
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo $FILESTEM2
	
	#Check whether a FASTQ file is properly formatted
    cutadapt -o /dev/null $file || error_exit "$LINENO: Fastq format error"
	
	#Remove Illumina Truseq adapters from all paired data using the 'regular 3'  paired adapters' function (these adapters are on the RHS of both R1 and R2 reads in all 8 samples)
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cores=12 -o ../1_trimmedData/$FILESTEM2"_R1_Iclean.fq.gz" -p ../1_trimmedData/$FILESTEM2"_R2_Iclean.fq.gz" $file $FILESTEM"_R2.fastq.gz" || error_exit "$LINENO: Error removing Illumina R1 adapters in $file"

	
	#Then remove Clontech SMARTer adapters from all paired data using 'linked adapters' function for paired end reads (No filtering for "anchored" adapters at this step - see notes, also not expected in G6 or G7)
	
	#as the length of target sequences is varied and unknown the possibility for the reverse completment of the second adapter is allowed by trimming as linked adapters.	
	#For linked adapters even though -a and -A are used (which indicate 3' adapters in the single line above) linked adapters assume the 5' is FIRST, and the 3' is SECOND.
	#Therefore to trim FORWARD R1 reads, they are called as <forward primer in 5'-3' orientation> ... <reverse complement of reverse primer> 
	#to trim REVERSE R2 reads, they are called as <reverse primer in 5'-3' orientation> ... <reverse complement of forward primer>	[https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters]
	#Reads shorter than 30bp are filtered out
	cutadapt -a AAGCAGTGGTATCAACGCAGAGTACNNNNN...AGTACTCTGCGTTGATACCACTGCTT -A AAGCAGTGGTATCAACGCAGAGTACT...NNNNNGTACTCTGCGTTGATACCACTGCTT --cores=12 --minimum-length 30 -o ../1_trimmedData/$FILESTEM2"_R1_clean.fq.gz" -p ../1_trimmedData/$FILESTEM2"_R2_clean.fq.gz" ../1_trimmedData/$FILESTEM2"_R1_Iclean.fq.gz" ../1_trimmedData/$FILESTEM2"_R2_Iclean.fq.gz" || error_exit "$LINENO: Error removing SMARTer adapters in $file"
done 
	

#run cleaned data through fastqc

