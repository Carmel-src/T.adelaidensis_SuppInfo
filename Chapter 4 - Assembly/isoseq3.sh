#!/bin/bash
#
#
# This script is used to characterise full length transcripts de novo using smrttools and bioconda as specified by PacBio workflow
# see: https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.0.md and follow links for bioconda references
#
# The paths in the script assume that a specific directory structure has been set up.
# This directory structure is for the Nectar cloud VM

# open a screen
# usage from /mnt
# bash isoseq3.sh <ctrl-a ctrl-d>
#
# Carmel Maher
# April 2019


#---------Installed Programs---------

# smrttools				-> smrttools-release_6.0.0.47835
	# help:
	# /pacbio/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/smrtcmds/bin/isoseq3 -h
	# /pacbio/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/smrtcmds/bin/isoseq3 --version
	# Version: isoseq3 3.0.0 (commit v3.0.0-7-gcc6cddd)

# python 2.7 			-> to run miniconda - automatically installed on nectar
# miniconda 			-> recommended python 2.7 compatible version by PacBio workflow. Miniconda2 | VER: 4.6.14
# ccs 					-> installed through conda, inside pbccs package | VER: pbccs-3.4.1
# lima					-> installed through conda | VER: lima-1.9.0


#------adjust these for your run------

#as per the local installation, the smrttools program directories are (in /mnt). Notable directories listed in installation: 
SMRT_ROOT=/mnt/pacbio/smrtlink

# Isoseq3 executeable location (from /mnt)
Iso3_DIR=$SMRT_ROOT/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/smrtcmds/bin

DATA=/mnt/data
MOVIE="MAH6260A1_m54196_190204_223227"
THREADS=0 	#0 for the Iseseq3 tools =autodetection of threads and is the default


#------------------------------------


function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}


# go to the working data directory
cd $DATA || error_exit "$LINENO: directory error"


# Generate consensus sequnces (ccs file) from raw subread data
# NOTE: this automatically runs with default autodetection of threads =  8
ccs $MOVIE.subreads.bam $MOVIE.ccs.bam --noPolish --minPasses 1 || error_exit "$LINENO: consensus sequence ccs error"

echo " - - - consensus sequnces generated"


# The primers.fasta as per the recommended Clontech SMARTer cDNA library prep - no barcodes used for this sample 
# primers.fasta has been uploaded to /mnt/data and does not need to be created here:
		# primers.fasta
		# >primer_5p
		# AAGCAGTGGTATCAACGCAGAGTACATGGG
		# >primer_3p
		# GTACTCTGCGTTGATACCACTGCTT
		
		
# Remove primers and demultiplex:
lima --isoseq --dump-clips --no-pbi $MOVIE.ccs.bam primers.fasta demux.bam || error_exit "$LINENO: Primer demultiplex lima error"

echo " - - - primers removed and demultiplexed"

			#!#!#!#!# 
      		# CHECK #  Note: A search using Git bash on the output files for any of the primers listed for Isoseq by Pacific biosciences returns no results after this script
      		#!#!#!#!#

# *********************
# From here on, execute the following steps for each output BAM file. -- only one in this case.
# isoseq3 tool usage: isoseq3 <tool>
# due to nectar permissions and installation location, usage: $Iso3_DIR/isoseq3 <tool>
# *********************


### cluster Tool: Cluster CCS reads and generate unpolished transcripts.
# recommended to give this step as many cores as possible
# Usage
# isoseq3 cluster [options] input output
# Example
# isoseq3 cluster <two types of potential inputfile> demux.bam <OR> $MOVIE.consensusreadset.xml  unpolished.bam

# Cluster consensus sequences to generate unpolished transcripts:
# note the demux.bam is named after the headers used in primers.fasta
$Iso3_DIR/isoseq3 cluster --verbose -j $THREADS demux.primer_5p--primer_3p.bam unpolished.bam || error_exit "$LINENO: Isoseq3 cluster error"

echo " - - - css reads clustered"


### polish Tool: Polish transcripts using subreads.
# Usage
# isoseq3 polish [options] input_1 input_2 output

# Polish transcripts using subreads:
$Iso3_DIR/isoseq3 polish --verbose -j $THREADS unpolished.bam $MOVIE.subreads.bam polished.bam || error_exit "$LINENO: Isoseq3 polish error"

echo " - - - transcripts polished"


### summarize Tool: Create a .csv-format barcode overview from transcripts.
# Usage
# isoseq3 summarize [options] input output
$Iso3_DIR/isoseq3 summarize --verbose polished.bam summary.csv || error_exit "$LINENO: Isoseq3 summary error"

echo " - - - summary file created"

echo " - - - done"