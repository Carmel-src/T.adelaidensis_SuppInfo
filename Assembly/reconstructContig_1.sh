#!/bin/bash
#
## this script is for coding region reconstruction of Cogent gene families as these must be completed individually.
# This script assumes you are following Running-Cogent.txt and using the same directory structure
#This script to be used after the process_kmer_to_graph.py step

# script in /mnt/IsoSeq-analysis/src
# usage from /mnt/IsoSeq-analysis/data/Cogent/hq
## bash ../../../src/reconstructContig.sh

# Carmel Maher
# August 2019

#Terry's error exit
function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message
    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

cd ~/mnt/IsoSeq-analysis/data/Cogent/hq

for d in */ ; do 

	DIR=${d%*}	
	echo $DIR
	
	reconstruct_contig.py -S T.adelaidensis ./$DIR || error_exit "$LINENO: Error "$DIR""
	
done
echo compete