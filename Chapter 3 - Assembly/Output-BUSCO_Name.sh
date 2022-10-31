#!/bin/bash
#
#
# These command lines were used to rename and copy assembled trinity & long read transcript datasets into the /BUSCO directory so all files have similar naming convention for input into the BUSCO script.
#
# The paths assume that a specific directory structure has been set up.
#
# Modules required: none
#
# 
# Carmel Maher
# August 2020 edited Oct 2022




#to grab the trinity output for BUSCO:
#-------------------------------------


###SHORT_READ ASSEMBLIES

# go to the working directory - working from /clean
# cd /scratch/user/mahe0050/DE-analysis/1_trimmedData/q  || error_exit "$LINENO: Directory Error 1"

#for file in *R1_cleanq.fq.gz
#	do
#
#		FILESTEM=${file%_*}
#		#this FILESTEM only cuts to _clean, the _R1 is included in stem
#		FILESTEM=${FILESTEM/R1/}
#		#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)

#All individual samples are named as below:
#		
#		/scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-sep/data/$FILESTEM"trinity-sep.Trinity".fasta

#Concatenated assembly of all samples is named:
#		
#		/scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2/trinity-all.Trinity.fasta /scratch/user/mahe0050/BUSCO/trinity-all_assembled.fasta  || error_exit "$LINENO: Error copying concatenated trinity-all output"
	
	
#####
#Note that using the above method all single sample trinity outputs can be called without editing or moving them.
#####


#-------------------------------------
### LONG READS (manually moved/renamed line by line)


cp -i /scratch/user/mahe0050/IsoSeq-analysis/data/Cogent/collected/hq.fasta.no5merge.collapsed.rep.fa /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep.fasta
	
		# -> one long-read file should end up in /scratch/user/mahe0050/IsoSeq-analysis/BUSCO named hq.fasta.no5merge.collapsed.rep_assembled.fasta This contains the full length of all non-redundant high quality transcripts

cp -i /scratch/user/mahe0050/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds /scratch/user/mahe0050/BUSCO/pygmy.ANGEL.cds.fasta
	
		# -> one long-read file should end up in /scratch/user/mahe0050/BUSCO named pygmy.ANGEL_assembled.fasta This contains the predicted coding sequence of non-redundant unique transcripts before protein clustering
		
cp -i /scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta /scratch/user/mahe0050/BUSCO/reference_transcripts.1Lv.clean.fasta
	
		# -> one long-read file should end up in /scratch/user/mahe0050/IsoSeq-analysis/BUSCO named reference_transcripts.1Lv.clean_assembled.fasta This contains the full length transcripts corresponding to the representative for clustered proteins translated from predicted open read frame


# header formats:
head /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep.fasta
>PB.1.1|000643|path0:1-1879(+)|transcript/14742 transcript/14742 full_length_coverage=3;length=1887;num_subreads=52

head /scratch/user/mahe0050/BUSCO/pygmy.ANGEL.cds.fasta
>PB.2.1|002537|path0:1-1624(+)|transcript/18304|m.1 type:likely-NA len:135 strand:+ pos:282-686

head /scratch/user/mahe0050/BUSCO/reference_transcripts.1Lv.clean.fasta
>PB.2.1|002537|path0:1-1624(+)|transcript/18304 transcript/18304 full_length_coverage=2;length=1627;num_subreads=26


grep -c ">" /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep.fasta
15729
grep -c ">" /scratch/user/mahe0050/BUSCO/pygmy.ANGEL.cds.fasta
13882
grep -c ">" /scratch/user/mahe0050/BUSCO/reference_transcripts.1Lv.clean.fasta
9813


# NOTE the files have symbols in faste headers such as "/" and "+" before transcript IDs which need to be removed for BUSCO

# after sed to remove all "/" resulted in an error at metaeuk headers "ValueError: could not convert string to float: '+'"
# after sed to remove all "+" same error

# Goal is a PERCENT score of alignments to BUSCOs and which specific transcripts match is no explored further. Therefore all headers were grossly simplified & truncated to remove symbols

cut -d '|' -f1 /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep.fasta > hq.fasta.no5merge.collapsed.rep_.fasta

cut -d '|' -f1 /scratch/user/mahe0050/BUSCO/pygmy.ANGEL.cds.fasta > pygmy.ANGEL.cds_.fasta

cut -d '|' -f1 /scratch/user/mahe0050/BUSCO/reference_transcripts.1Lv.clean.fasta > reference_transcripts.1Lv.clean_.fasta


# total number of sequences retained as counted above
grep -c ">" /scratch/user/mahe0050/BUSCO/hq.fasta.no5merge.collapsed.rep_.fasta
15729
grep -c ">" /scratch/user/mahe0050/BUSCO/pygmy.ANGEL.cds_.fasta
13882
grep -c ">" /scratch/user/mahe0050/BUSCO/reference_transcripts.1Lv.clean_.fasta
9813


# two of these files will have "duplicate" headers in this format as we know that some transcripts gave rise to >1 predicted cds. numbers were appended as below to make all sequence headers unique:
awk '/^>/{$0=$0"_"(++i)}1'  pygmy.ANGEL.cds_.fasta > pygmy.ANGEL.cds__.fasta
awk '/^>/{$0=$0"_"(++i)}1'  hq.fasta.no5merge.collapsed.rep_.fasta > hq.fasta.no5merge.collapsed.rep__.fasta

grep -c ">" hq.fasta.no5merge.collapsed.rep__.fasta
15729
grep -c ">" pygmy.ANGEL.cds__.fasta
13882
