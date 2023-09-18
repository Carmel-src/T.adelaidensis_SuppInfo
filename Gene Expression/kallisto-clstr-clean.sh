#!/bin/bash
#
#
# This script is used to estimate transcript counts using Kallisto - using HiSeq reads against an IsoSeq reference. 
# This reference was created from Iso-Seq data, clustered by CDHit (based on pygmy.ANGEL.pep, cutoff at .99, job 1598936109) and subset to include the representative transcript for each cluster using dplyR and seqtk.
# The paths in the script assume that a specific directory structure has been set up.
# Kallisto must be run on paired samples separately
# 
# usage on flinder's Deep Thought Machine via SLURM
# usage from within: conda activate ShortConda
#
# Carmel Maher & Terry Bertozzi
# Dec 2019, last altered Oct 2022


#------adjust these for your run-----

# module load Miniconda3/4.9.2
# conda activate ShortConda
# conda install kallisto
# kallisto, version 0.48.1


# ReferenceDIR=/scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk #contains reference_transcripts.1Lvclean.fasta
# KallistoEXE=/mnt/Prog/miniconda3/envs/anaKallisto/bin
KallistoDIR=/scratch/user/mahe0050/DE-analysis/2_alignedData/kallisto-clstr-cleanq
ReadDIR=/scratch/user/mahe0050/DE-analysis/1_trimmedData/q #contains trimmed R1, R2

#------------------------------------

function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

# *** make the index: ***
 cd /scratch/user/mahe0050/DE-analysis/2_alignedData/kallisto-clstr-cleanq

#create the index
#the below contains the representative transcript for clustered isoforms with (most) of the poly-a tails trimmed:
# kallisto index -i reference_transcripts.1Lv.clean.idx /scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta || error_exit "$LINENO: kallisto index error"

# ------------

## [build] loading fasta file /scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta
## [build] k-mer length: 31
## [build] counting k-mers ... done.
## [build] building target de Bruijn graph ...  done 
## [build] creating equivalence classes ...  done
## [build] target de Bruijn graph has 59186 contigs and contains 17770125 k-mers 

# ------------


# need to be in the input directory for for file
cd $ReadDIR

for file in *R1_cleanq.fq.gz
do
	FILESTEM=${file%_*}
	#this FILESTEM only cuts to _clean, the _R1 is included
	FILESTEM=${FILESTEM/R1/}
	#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)

echo $FILESTEM

kallisto quant -i $KallistoDIR/reference_transcripts.1Lv.clean.idx -o $KallistoDIR/$FILESTEM"kallisto-out" -b, --bootstrap-samples=100 --threads=12 --pseudobam $file $FILESTEM"R2_cleanq.fq.gz" || error_exit "$LINENO: kallisto error at $FILESTEM"


done

echo "kallisto-clstr-cleanq done"
