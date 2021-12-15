#!/bin/bash

# This script runs reconstruct_contig.py on every subdirectory contained in the top directory. reruns with a kmer parameter of +3 (to a max of 99) for failed directories.
# This script assumes directories and files exist as up until the Coding Genome Reconstruction step in the Running Cogent GitHub wiki 
# 
# usage: 	bash <path/script> <full directory containing hq_folders> <starting kmer value> <Species Name>
# e.g.	 	bash /mnt/IsoSeq-analysis/src/TBreconstructContig-edit.sh /mnt/IsoSeq-analysis/data/Cogent/hq 27 T.adelaidensis

# DIRNAME=$1
# species=$3

kmerSize=$2

find $1 -mindepth 1 -maxdepth 1 -type d | while read line; do #lists off all the /hq_* folders within the input directory and passes them through the script one by one (by line)
	FILE=$line/cogent2.fa
	if [ ! -f $FILE ]; then
		echo "failed no cogent2.fa Increase K-mer size" #Note all jobs will "fail" on first attempt unless you are running this in a previously reconstructed directory. In the latter case only the failed directories will continue to rerun.
		while [ ! -f "$FILE" ];do
			kmerSize=$((kmerSize + 3))
			rerunCMD="reconstruct_contig.py -k ${kmerSize} -S $3 $line" #reconstruction script increasing kmer parameter by 3 -note because all directories fail initially the input kmer value should be 3 below the desired actual value.
			echo ${rerunCMD}
			eval ${rerunCMD}
			if [ -f $FILE ]; then 
				touch $1/hq_kmerSize.txt #This text file will give a list of the used kmerSize per hq_* folder
				echo "SUCCESS! The increased Kmer size $kmerSize was successful for $line." >> $1/hq_kmerSize.txt
			fi
			if [ ${kmerSize} -gt 99 ];then
				cp $line/in.fa $line/cogent2.fa # if a folder reruns until 99kmer parameter length is reached the input is simply copied into cogent2.fa
				touch $1/failed-jobs.txt #This text file will remain empty or give a list of all failed folders
				echo "Failed to succeed reconstruction with largest Kmer on ${line}. Copied input transcripts to cogent2.fa as output." >> $1/failed-jobs.txt
			fi
		done
		kmerSize=$2
	fi
done

if [ -s $1/failed-jobs.txt ]; then
    echo "failed-jobs.txt not empty. Some jobs failed. Please re-run them."
else
    echo "failed-jobs.txt empty or doesnt exist. All jobs completed."
fi
