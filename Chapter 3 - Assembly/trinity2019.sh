#!/bin/bash
#
#
# This script is used to assemble transcripts using TRINITY
#
# The paths in the script assume that a specific directory structure has been set up.
#
# Modules required:
module add trinity/2.8.4
module add bowtie/2.3.5.1	
module add fastqc/0.11.8
module add OracleJava/12.0.1
module add samtools/1.9
#jellyfish
#salmon

# usage through slurm
#
# Carmel Maher & Terry Bertozzi
# Sep 2017 - edited Dec 2019


# how to incorporate singletons? https://www.biostars.org/p/149884/


#------adjust these for your run-----
THREADS=16

#Check available memory
# free -g
# max_memory
MEM="64G"

#------------------------------------


function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}


# go to the working directory - working from /clean?
cd /scratch/user/mahe0050/DE-analysis/1_trimmedData

# *** #

# To run trinity separately on each sample
mkdir -p ./trinity-sep
mkdir -p ./trinity-sep/temp

for file in *R1_clean.fq.gz
do
	FILESTEM=${file%_*}
	#this FILESTEM only cuts to _clean, the _R1 is included
	FILESTEM=${FILESTEM/R1/}
	#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)
	
	##concatenate relevant singleton file to R1 as library is not stranded?  
	cat ./$file ./$FILESTEM"singletons_clean.fq.gz" > ./trinity-sep/$FILESTEM"concat_R1.fq.gz" || error_exit "$LINENO: Error concatenating R1 and singletons"
	
	# assemble transcripts from a single paired sample file
	#note trinity requires output directory with "trinity" in name
	Trinity --seqType fq --max_memory $MEM --bflyCalculateCPU --left ./trinity-sep/$FILESTEM"concat_R1.fq.gz" --right ./$FILESTEM"R2_clean.fq.gz" --CPU $THREADS --output ./trinity-sep/$FILESTEM"trinity"/trinity || error_exit "$LINENO: Error running trinity-sep at $FILESTEM"

	
	#clean up
	rm -rf ./trinity-sep/$FILESTEM"concat_R1.fq.gz" || error_exit "$LINENO: Error removing trinity-sep concat at $FILESTEM"
	
done
echo "Trinity separate samples complete"

mkdir -p ./testes
mv G6T* ./testes
echo "testes moved out of working dir"

# *** #

#To create one master kidney transcript set, sample files to be concatenated and named appropriately
# create temp & working dir
mkdir -p ./trinity-all
mkdir -p ./trinity-all/temp

	#concatenate relevant sample files for one transcriptome - add singletons to R1 as its not stranded anyway? 
	cat ./*_R1_clean.fq.gz ./*_singletons_clean.fq.gz >> ./trinity-all/temp/PBTKI_cat_R1.fq.gz || error_exit "$LINENO: Error concatenating R1 and singletons"
	cat ./*_R2_clean.fq.gz >> ./trinity-all/temp/PBTKI_cat_R2.fq.gz || error_exit "$LINENO: Error concatenating R2"

	Trinity --seqType fq --max_memory $MEM --bflyCalculateCPU --left ./trinity-all/temp/PBTKI_cat_R1.fq.gz --right ./trinity-all/temp/PBTKI_cat_R2.fq.gz --CPU $THREADS --output ./trinity-all/trinity || error_exit "$LINENO: Error running Trinity-all"

	#clean up
	rm -rf ./trinity-all/temp || error_exit "$LINENO: Error removing trinity-all temp directory"
echo "Trinity concatenated kidney samples complete"