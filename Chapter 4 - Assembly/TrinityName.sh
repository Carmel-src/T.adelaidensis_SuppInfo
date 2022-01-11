#!/bin/bash
#
#
# This script is used to rename and copy trinity output files with their input file "filestem" into the /results directory
#
# The paths in the script assume that a specific directory structure has been set up.
#
# Modules required: none
#
# usage: 	<path to script>TrinityName.sh
# 			/home/users/cmaher/src/TrinityName.sh

# Carmel Maher
# October 2017


#------adjust these for your run-----

#------------------------------------


function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}


# go to the working directory - working from /clean
# cd /home/users/cmaher/data/clean

for file in *R1_clean.fq.gz
do

	FILESTEM=${file%_*}
	#this FILESTEM only cuts to _clean, the _R1 is included
	FILESTEM=${FILESTEM/R1/}
	#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)
	
	cp ../../results/$FILESTEM"trinity"/trinity/Trinity.fasta ../../results/$FILESTEM"trinity".fasta
	
done
	