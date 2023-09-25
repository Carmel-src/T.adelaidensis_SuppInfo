#!/bin/bash
# Please note that you need to adapt this script to your job
# Submitting as is will fail cause the job to fail 
# The keyword command for SLURM is #SBATCH --option
# Anything starting with a # is a comment and will be ignored
# ##SBATCH is a commented-out #SBATCH command
# SBATCH and sbatch are identical, SLURM is not case-sensitive
##################################################################
# Change FAN to your fan account name
# Change JOBNAME to what you want to call the job
# This is what is shows when attempting to Monitor / interrogate the job,
# So make sure it is something pertinent!
#
#SBATCH --job-name=mahe0050_BLAST-iso
#
##################################################################
# If you want email updates form SLURM for your job.
# Change MYEMAIL to your email address
#SBATCH --mail-user=carmel.maher@flinders.edu.au
#SBATCH --mail-type=ALL
# 
# Valid 'points of notification are': 
# BEGIN, END, FAIL, REQUEUE. 
# ALL means all of these
##################################################################
# Tell SLURM where to put the Job 'Output Log' text file. 
# This will aid you in debugging crashed or stalled jobs.
# You can capture both Standard Error and Standard Out
# %j will append the 'Job ID' from SLURM. 
# %x will append the 'Job Name' from SLURM 
# %
#SBATCH --output=/home/mahe0050/%x-%j.out.txt
#SBATCH --error=/home/mahe0050/%x-%j.err.txt
##################################################################
# The default partition is 'general'. 
# Valid partitions are general, gpu and melfu
##SBATCH --partition=PARTITIONNAME
#
##################################################################
# Tell SLURM how long your job should run for as a hard limit. 
# My setting a shorter time limit, it is more likely that your
# job will be scheduled when attempting to backfill jobs. 
# 
# The current cluster-wide limit is 14 Days from Start of Execution.
# The timer is only active while your job runs, so if you suspend
# or pause the job, it will stop the timer.
#
# The command format is as follows: #SBATCH --time=DAYS-HOURS
# There are many ways to specify time, see the SchedMD Slurm 
# manual pages for more. 
#SBATCH --time=0-12
#
##################################################################
# How many tasks is your job going to run? 
# Unless you are running something that is Parallel / Modular or
# pipelined, leave this as 1. Think of each task as a 'bucket of
# resources' that stand alone. Without MPI / IPC you can't talk to 
# another bucket!
#
#SBATCH --ntasks=1
#
# If each task will need more that a single CPU, then alter this 
# value. Remember, this is multiplicative, so if you ask for 
# 4 Tasks and 4 CPU's per Task, you will be allocated 16 CPU's 
#SBATCH --cpus-per-task=1
##################################################################
# Set the memory requirements for the job in MB. Your job will be
# allocated exclusive access to that amount of RAM. In the case it
# overuses that amount, Slurm will kill the job. The default value is 
# around 2GB per CPU you ask for.
#
# Note that the lower the requested memory, the higher the
# chances to get scheduled to 'fill in the gaps' between other
# jobs. Pick ONE of the below options. They are Mutually Exclusive.
# You can ask for X Amount of RAM per CPU (MB by default).
# Slurm understands K/M/G/T For Kilo/Mega/Giga/Tera Bytes.
#
#SBATCH --mem-per-cpu=4G
# Or, you can ask for a 'total amount of RAM'. If you have multiple 
# tasks and ask for a 'total amount' like below, then SLURM will 
# split the total amount to each task evenly for you.
##SBATCH --mem=G
##################################################################
# Change the number of GPU's required for you job. The most GPU's that can be 
# requested is 2 per node. As there are limited GPU slots, they are heavily 
# weighted against for Fairshare Score calculations. 
# You can request either a 'gpu:telsa_v100:X' or a 'gpu:x'
# 
# You can either request 0, or omit this line entirely if you 
# a GPU is not needed. 
#
##################################################################
# Load any modules that are required. This is exactly the same as 
# loading them manually, with a space-separated list, or you can 
# write multiple lines.
# You will need to uncomment these.

module load BLAST+/2.12 #BLAST+: 

##################################################################
# This example script assumes that you have already moved your 
# dataset to /scratch as part of your HPC Pre-Job preparations. 
# Its best to use the $TMP/$TMPDIR setup for you here
# to allow for the HPC to auto-clean anything you 
# leave behind by accident. 
# If you have a job-array and need a shared directory for 
# data on /local, you will need to manually cleanup that 
# directory as a part of your job script. 


##################################################################
# Enter the command-line arguments that you job needs to run. 

source ~/.bashrc 

#Script input:

# This script is used to run BLASTn to compare the concatenated TRINITY assembly (as the query) to the long-read Iso-seq transcripts (as a database) generated for T adelaidensis.
# The paths in the script assume that a specific directory structure has been set up.
#
# Modules required:
# module load BLAST+/2.12 #BLAST+: 

# usage through slurm
#
# Carmel Maher
# Aug 2023


#------------------------------------
function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}


#Create output folder
cd /scratch/user/mahe0050/DE-analysis/2_alignedData/
#mkdir ./BLAST || error_exit "$LINENO: Directory Error 0"
cd /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST

#reference to index
# /scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta

#assembled transcripts to query
# /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2/trinity-all.Trinity.fasta


# BLAST does not like '|' in FASTA database headers
# simplify the headers
sed 's/|/-/g' /scratch/user/mahe0050/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.clean.fasta > /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipe.fasta || error_exit "$LINENO: sed error1"
sed 's/path/p/g' /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipe.fasta > /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipep.fasta || error_exit "$LINENO: sed error2"
sed 's/transcript/t/g' /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipep.fasta > /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipept.fasta || error_exit "$LINENO: sed error3"

# Indexing the reference sequence
# makeblastdb -in [input database] -out [output database] -dbtype [database type] 

makeblastdb -in /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST/reference_transcripts.1Lv.c-nopipept.fasta -parse_seqids -dbtype nucl -out Reference_cluster_db || error_exit "$LINENO: Index error"


# BLASTn
# blastn -query [query file] -db [database file] -out [output file] -num_threads 4 -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 5

blastn -query /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2/trinity-all.Trinity.fasta -db Reference_cluster_db -out Trinity-all_IsoSeq_BLASTn_t5-h5 -num_threads 4 -outfmt 6 -evalue 1e-5 -max_target_seqs 5 -max_hsps 5 || error_exit "$LINENO: BLAST error"

blastn -query /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2/trinity-all.Trinity.fasta -db Reference_cluster_db -out Trinity-all_IsoSeq_BLASTn_t5-h5 -num_threads 4 -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 5 || error_exit "$LINENO: BLAST error"

blastn -query /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2/trinity-all.Trinity.fasta -db Reference_cluster_db -out Trinity-all_IsoSeq_BLASTn_t1-h1 -num_threads 4 -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 || error_exit "$LINENO: BLAST error"

echo "job complete"


# ----------
#Ran separately by commenting parts of code out:

#To sort BLASTn output for single best hit per query decided to use output from -evalue 1e-5 -max_target_seqs 5 -max_hsps 5 and sort manually
cd /scratch/user/mahe0050/DE-analysis/2_alignedData/BLAST
sort -k1,1 -k11,11g Trinity-all_IsoSeq_BLASTn_t5-h5 | sort --merge -u  -k1,1 -o Trinity-all_IsoSeq_BLASTn_t5-h5_Efilt.txt || error_exit "$LINENO: sort 1 error"
#This will sort by query sequence and p-value (assuming default format 6 output), then will select all first unique results (lowest p-value) for each query

#Then sort by target sequence, in order to easily pull out groups of transcripts that match the same reference for visualisation
sort -k2,2 Trinity-all_IsoSeq_BLASTn_t5-h5_Efilt.txt -o Trinity-all_IsoSeq_BLASTn_t5-h5_Efilt-Tsort.txt || error_exit "$LINENO: sort 2 error"

To determine how many target sequences received a hit:
cut -f 2 Trinity-all_IsoSeq_BLASTn_t5-h5 | sort | uniq | wc -l

# ----------

##################################################################
# Once you job has finished its processing, copy back your results 
# and ONLY the results to /scratch, then clean-up the temporary 
# working directory

# No need to cleanup $BGFS, SLURM handles the cleanup for you. 
# Just dont forget to copy out your results, or you will lose them!

##################################################################