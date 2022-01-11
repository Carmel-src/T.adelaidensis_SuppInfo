!/bin/bash

# November 2019
# Carmel Maher
# This script is to submit a BBduk script to cut adapters and quality on raw reads through Flinders' DeepThought HPC Slurm
# This script therefore assumes submission to Flinders Deepthought with an assumed directory structure and use of available modules


#------------------------------------
# Required Modules:
module add bbmap #/38.46
module add fastqc #/0.11.8

#adapters location :ref=/cm/shared/apps/bbmap/38.46/resources/adapters.fa

THREADS=8

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

#	in = /scratch/user/mahe0050/DE-analysis/0_rawData/
#	out = /scratch/user/mahe0050/DE-analysis/1_trimmedData

#------------------------------------

#change directories into raw sequencing directory 1
cd /scratch/user/mahe0050/DE-analysis/0_rawData/AGRF_CAGRF13871_HCGYYBCXY || error_exit "$LINENO: directory error 1"
mkdir ../../1_trimmedData/temp

for file in *R1.fastq.gz
do
	echo $file
	FILESTEM=${file%_*}
	echo $FILESTEM
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo $FILESTEM2
	
	#remove Illumina adapters from all paired data using bbduk
	bbduk.sh in=$file in2=$FILESTEM"_R2.fastq.gz" out=../../1_trimmedData/temp/$FILESTEM2"_R1_Iclean.fq.gz" out2=../../1_trimmedData/temp/$FILESTEM2"_R2_Iclean.fq.gz" outs=../../1_trimmedData/temp/$FILESTEM2"_singletons_Iclean.fq.gz" literal=GATCGGAAGAGCACA,AGATCGGAAGAGCGT ktrim=r k=15 mink=15 hdist=1 threads=$THREADS || error_exit "$LINENO: Error removing Illumina adapters from R1 or R2 - directory 1"
	
	#remove SMARTer PCR tags from all paired data using bbduk
	bbduk.sh in=../../1_trimmedData/temp/$FILESTEM2"_R1_Iclean.fq.gz" in2=../../1_trimmedData/temp/$FILESTEM2"_R2_Iclean.fq.gz" out=../../1_trimmedData/$FILESTEM2"_R1_clean.fq.gz" out2=../../1_trimmedData/$FILESTEM2"_R2_clean.fq.gz" outs=../../1_trimmedData/$FILESTEM2"_singletons_clean.fq.gz" literal=AAGCAGTGGTATCAACGCAGAGTACNNNNN,GTACTCTGCGTTGATACCACTGCTT ktrim=r k=15 mink=15 hdist=1 qtrim=r trimq=10 minlength=30 threads=$THREADS || error_exit "$LINENO: Error removing SMARTer adapters from R1 or R2 - directory 1"
	
	#remove SMARTer PCR tags from the singleton file using bbduk
	bbduk.sh in=../../1_trimmedData/temp/$FILESTEM2"_singletons_Iclean.fq.gz" out=../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" literal=AAGCAGTGGTATCAACGCAGAGTACNNNNN,GTACTCTGCGTTGATACCACTGCTT ktrim=r k=15 mink=15 hdist=1 qtrim=r trimq=20 minlength=30 threads=$THREADS || error_exit "$LINENO: Error removing SMARTer adapters from singletons - directory 1"
	
	#concatenate the singleton files 
	cat ../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" >> ../../1_trimmedData/$FILESTEM2"_singletons_clean.fq.gz" || error_exit "$LINENO: Error concatenating singletons - directory 1"
	
	
	#clean up: i.e. remove temporary directories
	rm ../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" || error_exit "$LINENO: Error deleting singleton intermediate - directory 1"
	
done 
	
	
#change directories into raw sequencing directory 2 and do it again for these samples (output to same trimmed directory for further analysis together)
cd /scratch/user/mahe0050/DE-analysis/0_rawData/AGRF_CAGRF20022_CDNJBANXX || error_exit "$LINENO: directory error 2"

for file in *R1.fastq.gz
do
	echo $file
	FILESTEM=${file%_*}
	echo $FILESTEM
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo $FILESTEM2
	
#(repeat of the above but edit directories & names: same output folder as all samples have unique names)
	#remove Illumina adapters from all paired data using bbduk
	bbduk.sh in=$file in2=$FILESTEM"_R2.fastq.gz" out=../../1_trimmedData/temp/$FILESTEM2"_R1_Iclean.fq.gz" out2=../../1_trimmedData/temp/$FILESTEM2"_R2_Iclean.fq.gz" outs=../../1_trimmedData/temp/$FILESTEM2"_singletons_Iclean.fq.gz" literal=GATCGGAAGAGCACA,AGATCGGAAGAGCGT ktrim=r k=15 mink=15 hdist=1 threads=$THREADS || error_exit "$LINENO: Error removing Illumina adapters from R1 or R2 - directory 2"
	
	#remove SMARTer PCR tags from all paired data using bbduk
	bbduk.sh in=../../1_trimmedData/temp/$FILESTEM2"_R1_Iclean.fq.gz" in2=../../1_trimmedData/temp/$FILESTEM2"_R2_Iclean.fq.gz" out=../../1_trimmedData/$FILESTEM2"_R1_clean.fq.gz" out2=../../1_trimmedData/$FILESTEM2"_R2_clean.fq.gz" outs=../../1_trimmedData/$FILESTEM2"_singletons_clean.fq.gz" literal=AAGCAGTGGTATCAACGCAGAGTACNNNNN,GTACTCTGCGTTGATACCACTGCTT ktrim=r k=15 mink=15 hdist=1 qtrim=r trimq=10 minlength=30 threads=$THREADS || error_exit "$LINENO: Error removing SMARTer adapters from R1 or R2 - directory 2"
	
	#remove SMARTer PCR tags from the singleton file using bbduk
	bbduk.sh in=../../1_trimmedData/temp/$FILESTEM2"_singletons_Iclean.fq.gz" out=../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" literal=AAGCAGTGGTATCAACGCAGAGTACNNNNN,GTACTCTGCGTTGATACCACTGCTT ktrim=r k=15 mink=15 hdist=1 qtrim=r trimq=20 minlength=30 threads=$THREADS || error_exit "$LINENO: Error removing SMARTer adapters from singletons - directory 2"
	
	#concatenate the singleton files 
	cat ../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" >> ../../1_trimmedData/$FILESTEM2"_singletons_clean.fq.gz" || error_exit "$LINENO: Error concatenating singletons - directory 2"
	
	
	#clean up: i.e. remove temporary directories
	rm ../../1_trimmedData/$FILESTEM2"_singletons1.fq.gz" || error_exit "$LINENO: Error deleting singleton intermediate - directory 2"
done 

cd /scratch/user/mahe0050/DE-analysis/1_trimmedData
rm -rf temp || error_exit "$LINENO: Error removing temp directory"


#run cleaned data through fastqc
mkdir -p /scratch/user/mahe0050/DE-analysis/1_trimmedData/FastQC-clean

for file in *R1_clean.fq.gz
do
	echo $file
	FILESTEM2=`echo $file | cut -d "_" -f1,2 --output-delimiter="_"`
	echo "FILESTEM2" $FILESTEM2
	
	fastqc --noextract --threads $THREADS -o ./FastQC-clean ./$FILESTEM2"_R1_clean.fq.gz" ./$FILESTEM2"_R2_clean.fq.gz" ./$FILESTEM2"_singletons_clean.fq.gz" || error_exit "$LINENO: Error with FastQC at "$FILESTEM2""
done



#primers check
#zcat G1_KI_HCGYYBCXY_ATCACG_L001_R1.fastq.gz | grep "AGATCGGAAGAGCAC" -m10
 #zcat G1_KI_HCGYYBCXY_ATCACG_L001_R1.fastq.gz | grep "GATCGGAAGAGCACA" -m10
 #zcat G1_KI_HCGYYBCXY_ATCACG_L001_R2.fastq.gz | grep "AGATCGGAAGAGCGT" -m10
