#!/bin/bash
#
#
# This script is used to estimate transcript counts using Kallisto - using HiSeq reads against an IsoSeq reference. 
# This reference was created from Iso-Seq data, clustered by CDHit (based on pygmy.ANGEL.pep, cutoff at .99, job 1598936109) and subset to include the representative transcript for each cluster using dplyR and seqtk.
# The paths in the script assume that a specific directory structure has been set up.
# Kallisto must be run on paired samples separately
# 
# usage on NECTAR
# usage from within: conda activate anaKallisto
# usage from /mnt/DE-analysis/bash
# activate screen
# bash kallisto-clstr-clean.sh &> /mnt/DE-analysis/2_alignedData/kallisto-clstr-clean/kallisto-clstr-clean-log.txt
#
# Carmel Maher & Terry Bertozzi
# Dec 2019, Altered Sep, Dec 2020. 12/1/21


#------adjust these for your run-----

# Miniconda3 4.7.12
#(changed PYTHONPATH from Cogent install -> also this is Python Python 3.7.4, that was Python 2.7)
# conda activate anaKallisto
# conda install -c bioconda kallisto
# kallisto, version 0.46.1


# ReferenceDIR=/mnt/IsoSeq-analysis/data/Seqtk #contains reference_transcripts.1Lvclean.fasta
# KallistoEXE=/mnt/Prog/miniconda3/envs/anaKallisto/bin
KallistoDIR=/mnt/DE-analysis/2_alignedData/kallisto-clstr-clean
ReadDIR=/mnt/DE-analysis/1_trimmedData #contains trimmed R1, R2, and singletons

#------------------------------------

# *** make the index: (done manually before script) ***

# cd /mnt/DE-analysis/2_alignedData
# mkdir kallisto-clstr-clean
# cd /mnt/DE-analysis/2_alignedData/kallisto-clstr-clean

#create the index
#the below contains the representative transcript for clustered isoforms with (most) of the poly-a tails trimmed:
# kallisto index -i reference_transcripts.1Lv.clean.idx /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta


# ______

# [build] loading fasta file /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta
# [build] k-mer length: 31
# [build] counting k-mers ... done.
# [build] building target de Bruijn graph ...  done
# [build] creating equivalence classes ...  done
# [build] target de Bruijn graph has 59186 contigs and contains 17770125 k-mers

# Worked better than kallisto-clstr-trim_test because this .idx did not need --make-unique to replace repeated names with unique names

# ______


#------------------------------------

function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

# need to be in the input directory for for file
cd $ReadDIR

for file in *R1_clean.fq.gz
do
	FILESTEM=${file%_*}
	#this FILESTEM only cuts to _clean, the _R1 is included
	FILESTEM=${FILESTEM/R1/}
	#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)

echo $FILESTEM

kallisto quant -i $KallistoDIR/reference_transcripts.1Lv.clean.idx -o $KallistoDIR/$FILESTEM"kallisto-out" -b, --bootstrap-samples=100 --threads=8 --pseudobam $file $FILESTEM"R2_clean.fq.gz" || error_exit "$LINENO: kallisto error at $FILESTEM"


done
echo "kallisto-clstr-clean done"
