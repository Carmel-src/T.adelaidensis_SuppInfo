#!/bin/bash
#
#
# This script is used to rename and copy trinity output files with their input file "filestem" into the /BUSCO directory so all files have similar naming convention for input into the BUSCO script.
#
# The paths in the script assume that a specific directory structure has been set up.
#
# Modules required: none
#
# usage: 	bash Output-BUSCO_Name2020.sh
# location:			cd /scratch/user/mahe0050/BUSCO/

# Carmel Maher
# August 2020


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
cd /scratch/user/mahe0050/DE-analysis/1_trimmedData  || error_exit "$LINENO: Directory Error 1"

for file in *R1_clean.fq.gz
	do

		FILESTEM=${file%_*}
		#this FILESTEM only cuts to _clean, the _R1 is included in stem
		FILESTEM=${FILESTEM/R1/}
		#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)
		
		cp -i /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-sep/$FILESTEM"trinity"/trinity/Trinity.fasta /scratch/user/mahe0050/BUSCO/$FILESTEM"t_assembled".fasta  || error_exit "$LINENO: Error copying separate trinity outputs"
		
			# -> all files should end up in /scratch/user/mahe0050/BUSCO named G#_XX.t_assembled.fasta
		
	done

	echo "Trinity-sep done"

cp -i /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/trinity/Trinity.fasta /scratch/user/mahe0050/BUSCO/trinity-all_assembled.fasta  || error_exit "$LINENO: Error copying concatenated trinity-all output"
	
		# -> one ANGEL file should end up in /scratch/user/mahe0050/BUSCO named trinity-all_assembled.fasta

	echo "Trinity-all done"
	
cp -i /scratch/user/mahe0050/IsoSeq-analysis/data/Cogent/collected/hq.fasta.no5merge.collapsed.rep.fa /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep_assembled.fasta  || error_exit "$LINENO: Error copying Cogent output"
	
		# -> one Cogent file should end up in /scratch/user/mahe0050/BUSCO named hq.fasta.no5merge.collapsed.rep_assembled.fasta

	echo "Cogent done"

cp -i /scratch/user/mahe0050/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds /scratch/user/mahe0050/BUSCO/pygmy.ANGEL_assembled.fasta  || error_exit "$LINENO: Error copying ANGEL output"
	
		# -> one ANGEL file should end up in /scratch/user/mahe0050/BUSCO named pygmy.ANGEL_assembled.fasta

	echo "ANGEL done"

	
	
# all file sizes checked after run = matches

# Trinity outputs removed from scratch 20/12/2020 to save space
