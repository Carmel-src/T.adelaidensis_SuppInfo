#!/bin/bash
#
#
# This script is used to assess transcriptome outputs against BUSCO lineages 
# The paths in the script assume that a specific directory structure has been set up.
# Analysis run on eRSA NECTAR machine
#
#Usage:
# 			conda activate conda-BUSCO
#			screen
# 			bash run-conda-BUSCO.bash 2>&1 | tee BUSCO-out.txt

# Carmel Maher
# February 2021


# See notes file associated with usage: BUSCO_all_v5.0.sh
# Inputs include:
# 	Short-read Trinity outputs for all 8 Kidney samples assembled from short reads
# 	Short-read Trinity output for one file of all 8 kidney sequencing concatenated and assembled into a single transcript file
# 	Isoseq3 long-read de novo assembly full list of hq transcripts
# 	Isoseq3 long-read de novo assembly list of non-redundant representative transcripts of peptide clusters, with additional poly-a tail trimming


#------adjust these for your run-----

LINEAGE="vertebrata_odb10"

#------------------------------------

function error_exit
{
    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

# (commented out because it failed after "done" - rerun with # )

#go to the directory containing trimmed files to pull ID names
cd /mnt/DE-analysis/1_trimmedData

for file in *R1_clean.fq.gz
	do
		FILESTEM=${file%_*}
		#this FILESTEM only cuts to _clean, the _R1 is included in stem
		FILESTEM=${FILESTEM/R1/}
		#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)
		
		mkdir /mnt/IsoSeq-analysis/BUSCO/$FILESTEM"BUSCOout" || error_exit "$LINENO: Error creating trinity-sep output directory at $FILESTEM"
		
		cd /mnt/IsoSeq-analysis/BUSCO/
		
		busco -i /mnt/DE-analysis/2_alignedData/trinity/trinity-sep/$FILESTEM"trinity"/trinity/Trinity.fasta -l $LINEAGE -f -c 12 -o $FILESTEM"BUSCOout" -m transcriptome || error_exit "$LINENO: Error running BUSCO at $FILESTEM"
		
		cd /mnt/DE-analysis/1_trimmedData
 	done

#

cd /mnt/IsoSeq-analysis/BUSCO/

mkdir ./trinity-all_BUSCOout  || error_exit "$LINENO: directory error at trinity-all"

busco -i /mnt/DE-analysis/2_alignedData/trinity/trinity-all/trinity/Trinity.fasta -l $LINEAGE -f -c 12 -o trinity-all_BUSCOout -m transcriptome || error_exit "$LINENO: Error running BUSCO at trinity-all"

#

cd /mnt/IsoSeq-analysis/BUSCO/

mkdir ./reference-transcripts_BUSCOout || error_exit "$LINENO: directory error at reference-transcripts"

busco -i reference_transcripts.1Lv.clean_.fasta -l $LINEAGE -f -c 12 -o reference-transcripts_BUSCOout -m transcriptome || error_exit "$LINENO: Error running BUSCO at reference-transcripts"

#

mkdir ./ANGEL.cds_BUSCOout || error_exit "$LINENO: directory error at reference-transcripts"

busco -i pygmy.ANGEL__.cds -l $LINEAGE -f -c 12 -o ANGEL.cds_BUSCOout -m transcriptome || error_exit "$LINENO: Error running BUSCO at ANGEL.cds"

#

mkdir ./hq-fasta_BUSCOout || error_exit "$LINENO: directory error at hq-fasta"

busco -i hq.fasta.no5merge.collapsed.rep.1L__.fa -l $LINEAGE -f -c 12 -o hq-fasta_BUSCOout -m transcriptome || error_exit "$LINENO: Error running BUSCO at hq-fasta"

#

echo "all BUSCOs complete"
