
#Counting the Transcripts assembled by Trinity:

cd /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-sep/data

grep -c ">" G1_KI_trinity-sep.Trinity.fasta
grep -c ">" G2_KI_trinity-sep.Trinity.fasta
grep -c ">" G3_KI_trinity-sep.Trinity.fasta
grep -c ">" G4_KI_trinity-sep.Trinity.fasta
grep -c ">" G5_KI_trinity-sep.Trinity.fasta
grep -c ">" G6_KI_trinity-sep.Trinity.fasta
grep -c ">" G7_KI_trinity-sep.Trinity.fasta
grep -c ">" G8_KI_trinity-sep.Trinity.fasta


# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G1_KI_trinity-sep.Trinity.fasta
# 48161
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G2_KI_trinity-sep.Trinity.fasta
# 45075
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G3_KI_trinity-sep.Trinity.fasta
# 37039
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G4_KI_trinity-sep.Trinity.fasta
# 40893
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G5_KI_trinity-sep.Trinity.fasta
# 38719
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G6_KI_trinity-sep.Trinity.fasta
# 185868
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G7_KI_trinity-sep.Trinity.fasta
# 229145
# (base) [mahe0050@hpc-head01 data]$ grep -c ">" G8_KI_trinity-sep.Trinity.fasta
# 43621



# *****



# Concatenated file of all 8.

cd /scratch/user/mahe0050/DE-analysis/2_alignedData/trinity-all/data2

grep -c ">" trinity-all.Trinity.fasta


# (base) [mahe0050@hpc-head01 data2]$ grep -c ">" trinity-all.Trinity.fasta
# 359197

