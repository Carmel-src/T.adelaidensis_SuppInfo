###runningANGEL


#Follow Isoseq3 and Running-Cogent.sh first

### ANGEL: Robust Open Reading Frame prediction ### https://github.com/PacificBiosciences/ANGEL

# This represents line by line instructions not an executeable script.
# This version has been simplified. Script parameters and file names are accurate but directories are an example as this data was run on a different HPC machine. However all results were then placed in ~IsoSeq-analysis/data/ANGEL

#The python dependencies for ANGEL are:

numpy
Biopython
scikit-learn
CD-HIT version 4.8.1
conda install -n anaCogent scikit-learn


## Dumb ORF prediction ##

# dumb_predict.py takes as input a FASTA file. It outputs the longest ORF in all frames - use to create a top training dataset.

# Usage:
dumb_predict.py <fasta_filename> <output_prefix> 
       [--min_aa_length MIN_AA_LENGTH]
       [--use_firstORF]
       [--use_rev_strand] [--cpus CPUS]
# By default, only the forward strand is used. This is especially true for PacBio transcriptome sequencing output that often already has correct strand. 

mkdir /mnt/IsoSeq-analysis/data/ANGEL
cd /mnt/IsoSeq-analysis/data/ANGEL


dumb_predict.py /mnt/IsoSeq-analysis/data/Cogent/collected/hq.fasta.no5merge.collapsed.rep.fa pygmy.dumb --min_aa_length 100 --cpus 24

#------
13335  finished       7611  clusters

#------


## Creating a non-redundant training dataset ##
## ANGEL classifier training ##


angel_make_training_set.py pygmy.dumb.final pygmy.dumb.final.training --random --cpus 24


angel_train.py pygmy.dumb.final.training.cds pygmy.dumb.final.training.utr pygmy.dumb.final.classifier.pickle --cpus 12


## Robust ORF prediction ##
# based on both the ANGEL Classifier training and the dumb ORF prediction

# Sequences are output as fasta files with  whether they are complete, 5' partial, 3' partial, or internal in the header information, which also includes the aa length. 
# Preidctions are tagges as confident, likely, or suspicious, and dumb ORF predictions as dumb.
# The proportion of each isoform that is untranslated, as well as the position of the cds sequence is also output.

angel_predict.py ../Cogent/collected/hq.fasta.no5merge.collapsed.rep.fa pygmy.dumb.final.classifier.pickle pygmy --output_mode=best --min_angel_aa_length 100 --min_dumb_aa_length 100 --cpus 48

Output is written to pygmy.ANGEL.cds, pygmy.ANGEL.pep, pygmy.ANGEL.utr


#Note: The reverse strand options within ANGEL are not necessary on Isoseq data
