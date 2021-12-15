# This document contains notes and single lines run to run BUSCO 5.0 installed using the conda package on all short read and long read assembled transcript files.

# Conda Version 4.8.3 has been used for all analyses up to this point.

#update conda as version 4.8.4 or higher is recommended.
conda update -n base conda

#This analysis was conducted on eRSA's NECTAR research cloud


# Working directory of full (single line) hq transcript Isoseq dataset, and trimmed, non-redundant file of transcripts identified as cluster representatives.
/mnt/IsoSeq-analysis/data/Seqtk/hq.fasta.no5merge.collapsed.rep.1L.fa
/mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta

# Working directory of assembled (using TRINITY) short read transcript files
/mnt/DE-analysis/2_alignedData/trinity

	# concatenated short reads all assembled using trinity into largest transcript file/mnt/DE-analysis/2_alignedData/trinity/trinity-all/trinity
	#	/mnt/DE-analysis/2_alignedData/trinity/trinity-all/trinity

	# Then under 
	#	/trinity-sep/<sample_ID_trinity>/trinity

	# Sample IDs folders:
	#	G1_KI_trinity
	#	G2_KI_trinity
	#	G3_KI_trinity
	#	G4_KI_trinity
	#	G5_KI_trinity
	#	G6K_CDNJBANXX_trinity
	#	G7K_CDNJBANXX_trinity
	#	G8_KI_trinity

# -----
# from within /mnt/DE-analysis/1_trimmedData
#the following command will run through these ID names:

		# for file in *R1_clean.fq.gz
		# 	do

		# 		FILESTEM=${file%_*}
				#this FILESTEM only cuts to _clean, the _R1 is included in stem
		# 		FILESTEM=${FILESTEM/R1/}
				#removes R1 from FILESTEM (FILESTEM ends in _ therefore not needed in "text" names)
				
		# 		command for: /mnt/DE-analysis/2_alignedData/trinity/trinity-sep/$FILESTEM"trinity"

		# 	done

# -----

# *** ***
# Sample format compatibility

cd /mnt/IsoSeq-analysis/BUSCO/

# NOTE the Isoseq files have "/" in fasta headers before transcript IDs which need to be removed for BUSCO

# header formats:
head /mnt/IsoSeq-analysis/data/Seqtk/hq.fasta.no5merge.collapsed.rep.1L.fa
>PB.4.2|004079|path1:6-612(+)|transcript/24311 transcript/24311 full_length_coverage=2;length=640;num_subreads=60

head /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta
>PB.11.2|004815|path2:1047-3035(+)|transcript/13401 transcript/13401 full_length_coverage=2;length=1995;num_subreads=8

head /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds
>PB.11.2|004815|path2:1047-3035(+)|transcript/13401|m.5 type:dumb-complete len:311 strand:+ pos:214-1146

grep -c ">" /mnt/IsoSeq-analysis/data/Seqtk/hq.fasta.no5merge.collapsed.rep.1L.fa
15729
grep -c ">" /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta
9813
grep -c ">" /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds
13882

# ***
# after sed to remove all "/" resulted in an error at metaeuk headers "ValueError: could not convert string to float: '+'"
# after sed to remove all "+" same error

# Goal is a %score of alignments to BUSCOs and which transcripts match is no explored further. therefore all headers were grossly simplified

cut -d '|' -f1 /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta > reference_transcripts.1Lv.clean_.fasta

cut -d '|' -f1 /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds > pygmy.ANGEL_.cds

cut -d '|' -f1 /mnt/IsoSeq-analysis/data/Seqtk/hq.fasta.no5merge.collapsed.rep.1L.fa > hq.fasta.no5merge.collapsed.rep.1L_.fa

# total number of sequences retained as counted above
grep -c ">" reference_transcripts.1Lv.clean_.fasta
9813
grep -c ">" pygmy.ANGEL_.cds
13882
grep -c ">" hq.fasta.no5merge.collapsed.rep.1L_.fa
15729

# two of these files will have "duplicate" headers in this format. numbers appended as below:
awk '/^>/{$0=$0"_"(++i)}1'  pygmy.ANGEL_.cds > pygmy.ANGEL__.cds
awk '/^>/{$0=$0"_"(++i)}1'  hq.fasta.no5merge.collapsed.rep.1L_.fa > hq.fasta.no5merge.collapsed.rep.1L__.fa

grep -c ">" pygmy.ANGEL__.cds
13882
grep -c ">" hq.fasta.no5merge.collapsed.rep.1L__.fa
15729


# *** *** ***
# BUSCO Installation & usage

cd /mnt/IsoSeq-analysis/BUSCO

#Create a conda environment for the BUSCO Install

conda create -n conda-BUSCO -c bioconda -c conda-forge busco=5.0.0
conda activate conda-BUSCO

#The basic repeated usage of BUSCO for all samples listed below will be:
busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m transcriptome


# note: E-value cutoff for BLAST searches Default: 1e-03

# busco --list-datasets
#check $LINEAGE in script file
#		- vertebrata_odb10
#		may also --auto-lineage-euk

# script currently set to threads (-c) = 12

# *** ***
#Run BUSCO -> separate script run-conda-BUSCO.bash

# location of usage doesnt matter as paths are specified within script.
# use from within /mnt/IsoSeq-analysis/BUSCO/ in case of log files.
# Usage:
 			conda activate conda-BUSCO
			screen
 			bash run-conda-BUSCO.bash 2>&1 | tee BUSCO-out.txt


# *** ***
# generate plot.py and r script of BUSCO summaries

mkdir BUSCO_summaries

cp G1_KI_BUSCOout/short_summary.*_odb10.G1_KI_BUSCOout.txt BUSCO_summaries/.
cp G2_KI_BUSCOout/short_summary.*_odb10.G2_KI_BUSCOout.txt BUSCO_summaries/.
cp G3_KI_BUSCOout/short_summary.*_odb10.G3_KI_BUSCOout.txt BUSCO_summaries/.
cp G4_KI_BUSCOout/short_summary.*_odb10.G4_KI_BUSCOout.txt BUSCO_summaries/.
cp G5_KI_BUSCOout/short_summary.*_odb10.G5_KI_BUSCOout.txt BUSCO_summaries/.
cp G6K_CDNJBANXX_BUSCOout/short_summary.*_odb10.G6K_CDNJBANXX_BUSCOout.txt BUSCO_summaries/.
cp G7K_CDNJBANXX_BUSCOout/short_summary.*_odb10.G7K_CDNJBANXX_BUSCOout.txt BUSCO_summaries/.
cp G8_KI_BUSCOout/short_summary.*_odb10.G8_KI_BUSCOout.txt BUSCO_summaries/.
cp reference-transcripts_BUSCOout/short_summary.*_odb10.reference-transcripts_BUSCOout.txt BUSCO_summaries/.
cp ANGEL.cds_BUSCOout/short_summary.*_odb10.ANGEL.cds_BUSCOout.txt BUSCO_summaries/.
cp hq-fasta_BUSCOout/short_summary.*_odb10.hq-fasta_BUSCOout.txt BUSCO_summaries/.
cp trinity-all_BUSCOout/short_summary.*_odb10.trinity-all_BUSCOout.txt BUSCO_summaries/.

# (The path below is where the config files for conda-BUSCO environment are located)

python3 /mnt/Prog/miniconda3/envs/conda-BUSCO/bin/generate_plot.py -wd /mnt/IsoSeq-analysis/BUSCO/BUSCO_summaries

# *** *** *** 

mkdir BUSCO_SR_summaries

cp G1_KI_BUSCOout/short_summary.*_odb10.G1_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp G2_KI_BUSCOout/short_summary.*_odb10.G2_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp G3_KI_BUSCOout/short_summary.*_odb10.G3_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp G4_KI_BUSCOout/short_summary.*_odb10.G4_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp G5_KI_BUSCOout/short_summary.*_odb10.G5_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp G6K_CDNJBANXX_BUSCOout/short_summary.*_odb10.G6K_CDNJBANXX_SR_BUSCOout.txt BUSCO_summaries/.
cp G7K_CDNJBANXX_BUSCOout/short_summary.*_odb10.G7K_CDNJBANXX_SR_BUSCOout.txt BUSCO_summaries/.
cp G8_KI_BUSCOout/short_summary.*_odb10.G8_KI_BUSCOout.txt BUSCO_SR_summaries/.
cp trinity-all_BUSCOout/short_summary.*_odb10.trinity-all_BUSCOout.txt BUSCO_SR_summaries/.

python3 /mnt/Prog/miniconda3/envs/conda-BUSCO/bin/generate_plot.py -wd /mnt/IsoSeq-analysis/BUSCO/BUSCO_SR_summaries

# *** *** *** 

mkdir BUSCO_LR_summaries

cp reference-transcripts_BUSCOout/short_summary.*_odb10.reference-transcripts_BUSCOout.txt BUSCO_LR_summaries/.
cp ANGEL.cds_BUSCOout/short_summary.*_odb10.ANGEL.cds_BUSCOout.txt BUSCO_LR_summaries/.
cp hq-fasta_BUSCOout/short_summary.*_odb10.hq-fasta_BUSCOout.txt BUSCO_LR_summaries/.

python3 /mnt/Prog/miniconda3/envs/conda-BUSCO/bin/generate_plot.py -wd /mnt/IsoSeq-analysis/BUSCO/BUSCO_LR_summaries