# Running Cogent
# 2019

# Some commentary based on noted from Tessa Bradford & Terry Bertozzi, and the Cogent GitHub Tutorial page
# https://github.com/Magdoll/Cogent/wiki/Running-Cogent


# *****


# activate screen
screen -s cogent

# activate the environment
conda activate anaCogent

# make Cogent output directory
mkdir /mnt/IsoSeq-analysis/data/Cogent

# go to analysis working directory
cd /mnt/IsoSeq-analysis/data/Cogent

# from analysis directory the input data directory is ../
# ../2019-04-30_out/polished.hq.fasta

### the first output will be in the original data directory

# check paths (During installation, all paths were put into .profile)
echo $PATH
echo $PYTHONPATH
echo $LD_LIBRARY_PATH


# *****


#### Running Cogent ####
# https://github.com/Magdoll/Cogent/wiki/Running-Cogent#findsmall


### Family Finding ###
## Running Family Finding for a small dataset ##

# create k-mer profile of the input and calculate pairwise distances
# note: default K-mer sixe =30
python /mnt/Prog/Cogent/Cogent/run_mash.py --cpus=7 /mnt/IsoSeq-analysis/data/2019-04-30_out/polished.hq.fasta

# the output will be <fasta_filename>.k<sketch_size>.dist
## Output written to: /mnt/IsoSeq-analysis/data/2019-04-30_out/polished.hq.fasta.s1000k30.dist

# process the distance file and create the partitions for each gene family
process_kmer_to_graph.py ../2019-04-30_out/polished.hq.fasta ../2019-04-30_out/polished.hq.fasta.s1000k30.dist ./hq hq

# This would generate an output log file hq.partition.txt and for each partition ("gene" family or isoform set), a subdirectory called hq/<partition_number> which contains the subset of fasta seqeuences belonging to that family. Note that sequences that don't belong to any partition (because it has no similarity with other sequences) will be "unassigned" and noted in the partition log file


### Coding Genome Reconstruction ###

## each "gene" family must be done individually for coding region ##
# Reconstructed contigs will contain the whole coding region. Reconstructed file is called cogent2.renamed.fasta

# The script for an individual folder is:
# reconstruct_contig.py -S T.adelaidensis hq/hq_0

# Some testing required to include 5 folders which failed to run with kmer length 30 https://github.com/Magdoll/Cogent/issues/40

# modified script with a loop for failed folders based on the above github issue. #the loop is awkward so inputting 27 will start all at kmer size 30.
bash /mnt/IsoSeq-analysis/src/TBreconstructContig-edit.sh /mnt/IsoSeq-analysis/data/Cogent/hq 27 T.adelaidensis


# number of completed gene families
ls ./hq/*/COGENT.DONE | wc -l
4459
# OR
ls ./hq/*/cogent2.renamed.fasta | wc -l
4459

### Finding failed reconstruction jobs ###
ls hq | wc -l #number of items in ./hq
4460
# the added 1 is the text file reporting the successful kmer size for each "gene" family folder


# The following will list the names of folders that have not finished (Terry):
comm -3 <(find /mnt/IsoSeq-analysis/data/Cogent/hq -iname 'COGENT.DONE' -printf '%h\n' | sort -u) <(find /mnt/IsoSeq-analysis/data/Cogent/hq -maxdepth 1 -mindepth 1 -type d | sort) | sed -e 's/^.*hq\///'


#------
#see failed-jobs.txt and hq_kmerSize.txt
#------


# *****



#### Tutorial: Using Cogent to collapse redundant transcripts in absence of genome ####
# https://github.com/Magdoll/Cogent/wiki/Tutorial:-Using-Cogent-to-collapse-redundant-transcripts-in-absence-of-genome

### Creating the "fake genome" ###
## list of and number of unassigned sequences ##

# make sure you are in the correct Cogent directory and can see the file hq.partition.txt
tail -n 1 hq.partition.txt |tr ',' '\n' > unassigned.list
# creates the unassigned.list file in the same directory

tail -n 1 hq.partition.txt |tr ',' '\n' | wc -l
3193    #nb this number did not change with the addition of the 5 failed reconstruction folders. i.e. they failed completely and did not end up in the unassigned dataset. Now that they are fixed this dataset is unaffected.

# Make the unassigned hq sequences into a fasta file
# Make sure you are in the anaCogent environment so Cupcake works
# make sure the .py directory for cupcake is in the PATH
export PATH=$PATH:/mnt/Prog/cDNA_Cupcake/sequence
get_seqs_from_list.py ../2019-04-30_out/polished.hq.fasta unassigned.list > unassigned.fasta


##Concatenate unassigned with Cogent contigs ##

# put the reconstructed genes plus the unassigned single hq isoforms into a single fasta file
mkdir collected
cd collected
cat ../hq/*/cogent2.renamed.fasta ../unassigned.fasta > cogent_fake_genome.fasta


## Collapsing redundant isoforms ##
# create a SAM alignment. This can be done either with Minimap2 or GMAP
# We have used minimap2/
# Obtain a final set of unique (non-redundant) transcript Isoforms that can be used as a reference gene set
# as there is natural 5' degradation in RNA some sequences will represent identical isoforms which may not all be identified in clustering

# create an aligned, sorted SAM file
# these parameters are default or suggested in the Cogent tutorial page
export PATH=$PATH:/mnt/Prog/minimap2
minimap2 -ax splice -t 30  -uf --secondary=no cogent_fake_genome.fasta ../../2019-04-30_out/polished.hq.fasta > hq.fasta.sam

#------
[M::mm_idx_gen::2.392*1.00] collected minimizers
[M::mm_idx_gen::3.526*1.32] sorted minimizers
[M::main::3.528*1.32] loaded/built the index for 9960 target sequence(s)
[M::mm_mapopt_update::3.708*1.30] mid_occ = 31
[M::mm_idx_stat] kmer size: 15; skip: 5; is_hpc: 0; #seq: 9960
[M::mm_idx_stat::3.810*1.29] distinct minimizers: 6267388 (82.12% are singletons); average occurrences: 1.346; average spacing: 2.846
[M::worker_pipeline::13.901*5.65] mapped 25117 sequences
[M::main] Version: 2.11-r797
[M::main] CMD: /mnt/Prog/pacbio/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/private/thirdparty/minimap2/minimap2_2.11/binwrap/../../../../../private/thirdparty/minimap2/minimap2_2.11/bin/minimap2 -ax splice -t 30 -uf --secondary=no cogent_fake_genome.fasta ../../2019-04-30_out/polished.hq.fasta
[M::main] Real time: 13.942 sec; CPU: 78.558 sec
#------


# Then follow the collapse tutorial from Cupcake
# https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step
# there is another cupcake page with information but the example scripts are provided on the cogent tutorial page

# sort the SAM file
sort -k 3,3 -k 4,4n hq.fasta.sam > hq.fasta.sorted.sam

# Colapse identical isoforms to obtain a list of full length, unique, hq isoforms to use as reference transcripts

cd /mnt/Prog/cDNA_Cupcake/sequence
which collapse_isoforms_by_sam.py
/home/ubuntu/miniconda2/envs/anaCogent/bin/collapse_isoforms_by_sam.py


cd /mnt/IsoSeq-analysis/data/Cogent/collected

# collapse usage is:
# usage: collapse_isoforms_by_sam.py [-h] [--input INPUT] [--fq] -s SAM -o
#                                    PREFIX [-c MIN_ALN_COVERAGE]
#                                    [-i MIN_ALN_IDENTITY]
#                                    [--max_fuzzy_junction MAX_FUZZY_JUNCTION]
#                                    [--flnc_coverage FLNC_COVERAGE]
#                                    [--dun-merge-5-shorter]

collapse_isoforms_by_sam.py --input ../../2019-04-30_out/polished.hq.fasta -s hq.fasta.sorted.sam -c 0.94 -i 0.85 --dun-merge-5-shorter -o hq.fasta.no5merge
# these parameters taken from pers. comm Tessa Bradford, SA Museum.
# output is the name 'stem' the collapsed.rep.fa is for annotation
# The output files are <-o>.collapsed.gff, <-o>.collapsed.rep.fq, and <-o>.collapsed.group.txt.
# The naming system for the post-collapse isoform is PB.<loci_index>.<isoform_index>


#------
# more terminal output above excluded + example:

Ignored IDs written to: hq.fasta.no5merge.ignored_ids.txt
Output written to:
hq.fasta.no5merge.collapsed.gff
hq.fasta.no5merge.collapsed.group.txt
hq.fasta.no5merge.collapsed.rep.fa
Namespace(allow_extra_5exon=False, flnc_coverage=-1, fq=False, input='../../2019-04-30_out/polished.hq.fasta', max_3_diff=100, max_5_diff=1000, max_fuzzy_junction=5, min_aln_coverage=0.94, min_aln_identity=0.85, prefix='hq.fasta.no5merge', sam='hq.fasta.sorted.sam')

ls

cogent_fake_genome.fasta                       hq.fasta.no5merge.collapsed.rep.fa
hq.fasta.no5merge.collapsed.gff                hq.fasta.no5merge.ignored_ids.txt
hq.fasta.no5merge.collapsed.gff.unfuzzy        hq.fasta.sam
hq.fasta.no5merge.collapsed.group.txt          hq.fasta.sorted.sam
hq.fasta.no5merge.collapsed.group.txt.unfuzzy
#------


# hq.fasta.no5merge.collapsed.group.txt names the isoforms as PB.<loci_index>.<isoform_index> and lists the sollapsed identical isoforms
# Each locus consists of a strand-specific locus with isoforms that overlap by at least 1 bp. So PB.11.1 and PB.11.2 means this locus has two isoforms.

wc -l hq.fasta.no5merge.collapsed.group.txt
15729 hq.fasta.no5merge.collapsed.group.txt # loci (transcripts) found

get_abundance_post_collapse.py hq.fasta.no5merge.collapsed /mnt/IsoSeq-analysis/data/2019-04-30_out/polished.cluster_report.csv
WARNING: isoseq3 format detected. Output `length` column will be `NA`.
#NOTE this is true: hq.fasta.no5merge.collapsed.read_stat.txt has no length information and a whole column of NA values, however the headers of hq.fasta.no5merge.collapsed.rep.fa do still have length information
Read stat file written to hq.fasta.no5merge.collapsed.read_stat.txt
Abundance file written to hq.fasta.no5merge.collapsed.abundance.txt



# Get count information from the abundance and group .txt files
echo -e "pbid\tcount_fl" > output.collapsed.abundance.txt
#adds the headers

paste <(awk '{print $1}' hq.fasta.no5merge.collapsed.group.txt) <(awk -F , '{print NF}' hq.fasta.no5merge.collapsed.group.txt) >> output.collapsed.abundance.txt
#adds the PB.<loci_index>.<isoform_index> and number of isoforms



#longest transcript for each set of isoforms is in 
hq.fasta.no5merge.collapsed.rep.fa


### Move on to Running-ANGEL Instructions ###

