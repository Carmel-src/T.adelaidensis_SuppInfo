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
#SBATCH --job-name=mahe0050_Trinity
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
#SBATCH --time=14-0
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
#SBATCH --cpus-per-task=16
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
##SBATCH --mem-per-cpu=12G
# Or, you can ask for a 'total amount of RAM'. If you have multiple 
# tasks and ask for a 'total amount' like below, then SLURM will 
# split the total amount to each task evenly for you.
#SBATCH --mem=128G
##################################################################
# Change the number of GPU's required for you job. The most GPU's that can be 
# requested is 2 per node. As there are limited GPU slots, they are heavily 
# weighted against for Fairshare Score calculations. 
# You can request either a 'gpu:telsa_v100:X' or a 'gpu:x'
# 
# You can either request 0, or omit this line entirely if you 
# a GPU is not needed. 
#
#SBATCH --gres="gpu:0"
##################################################################
# Load any modules that are required. This is exactly the same as 
# loading them manually, with a space-separated list, or you can 
# write multiple lines.
# You will need to uncomment these.
 
module load singularity/3.6.3

##################################################################
# This example script assumes that you have already moved your 
# dataset to /scratch as part of your HPC Pre-Job preparations. 
# Its best to use the $TMP/$TMPDIR setup for you here
# to allow for the HPC to auto-clean anything you 
# leave behind by accident. 
# If you have a job-array and need a shared directory for 
# data on /local, you will need to manually cleanup that 
# directory as a part of your job script. 

# Example using the SLURM $BGFS Variable (the Parallel Filesystem)

cd $BGFS
mkdir $BGFS/data/

cd /scratch/user/mahe0050/DE-analysis/1_trimmedData/q

##################################################################
# Enter the command-line arguments that you job needs to run. 


# Carmel Maher & Terry Bertozzi
# Sep 2017 - last edited oct 2022


function error_exit
{
    # Exit function due to fatal error
    # Accepts 1 arg:
    # string - descriptive error message

    echo "${PROGNAME}: ${1:-"Unknown error"}" 1>&2
    exit 1
}

#Slurm script places you in a temporary directory, this shall be treated as the same 'level' as /scratch/user/mahe0050/DE-analysis/

# go to the working directory - working from $BGFS/data
pwd

# *** #

# To run trinity separately on each sample
	
	# assemble transcripts from a single paired sample file
	#note trinity requires output directory with "trinity" in name
	
# ***G1***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G1_KI_R1_cleanq.fq.gz  \
          --right /data/G1_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G1_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G1_KI_"

echo "Trinity G1 complete"

	
# ***G2***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G2_KI_R1_cleanq.fq.gz  \
          --right /data/G2_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G2_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G2_KI_"

echo "Trinity G2 complete"

	
# ***G3***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G3_KI_R1_cleanq.fq.gz  \
          --right /data/G3_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G3_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G3_KI_"

echo "Trinity G3 complete"

	
# ***G4***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G4_KI_R1_cleanq.fq.gz  \
          --right /data/G4_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G4_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G4_KI_"

echo "Trinity G4 complete"

	
# ***G5***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G5_KI_R1_cleanq.fq.gz  \
          --right /data/G5_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G5_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G5_KI_"

echo "Trinity G5 complete"

	
# ***G6***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G6_KI_R1_cleanq.fq.gz  \
          --right /data/G6_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G6_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G6_KI_"

echo "Trinity G6 complete"

	
# ***G7***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G7_KI_R1_cleanq.fq.gz  \
          --right /data/G7_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G7_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G7_KI_"

echo "Trinity G7 complete"

	
# ***G8***	

	singularity exec --home $BGFS/data:/home -B /scratch/user/mahe0050/DE-analysis/1_trimmedData/q:/data -e /scratch/user/mahe0050/DE-analysis/trinityrnaseq.v2.14.0.simg  Trinity \
          --seqType fq \
          --left /data/G8_KI_R1_cleanq.fq.gz  \
          --right /data/G8_KI_R2_cleanq.fq.gz \
          --verbose --CPU 16 --max_memory 128G \
          --output G8_KI_trinity-sep || error_exit "$LINENO: Error running trinity-sep at G8_KI_"

echo "Trinity G8 complete"


##############
# REMOVE all * , relative file positions, FILESTEM & other variables, they are not read within the container.
##############


##################################################################
# Once you job has finished its processing, copy back your results 
# and ONLY the results to /scratch, then clean-up the temporary 
# working directory
# This command assumes that the destination exists

mkdir /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-sep

cp -r /$BGFS/data /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-sep

# No need to cleanup $BGFS, SLURM handles the cleanup for you. 
# Just dont forget to copy out your results, or you will lose them!

##################################################################