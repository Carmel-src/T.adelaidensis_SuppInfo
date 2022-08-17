
#Counting the Transcripts assembled by Trinity:

# First 6 - including only the first sequencing run and that concatenated file of 6 samples used for troubleshooting.
cd /e/@DATA/@@PBT_KidneyHiSeq_WorkingData_C.Maher/Not used or redundant/2019 Hiseq folders TANGO Trinity-Corset pipe/results/trinity.fasta

$ grep -c ">" G1_KI_HCGYYBCXY_ATCACG_L001_trinity.fasta
#47086

$ grep -c ">" G2_KI_HCGYYBCXY_CGATGT_L001_trinity.fasta
#43621

$ grep -c ">" G3_KI_HCGYYBCXY_TTAGGC_L001_trinity.fasta
#35960

$ grep -c ">" G4_KI_HCGYYBCXY_TGACCA_L001_trinity.fasta
#39482

$ grep -c ">" G5_KI_HCGYYBCXY_ACAGTG_L001_trinity.fasta
#37474

$ grep -c ">" G8_KI_HCGYYBCXY_GCCAAT_L001_trinity.fasta
#43456

$ grep -c ">" concat_trinity.fasta
#88654




# *****



# All 8 Samples main TRINITY run - including concatenated file of all 8.

pwd
/e/@DATA/@@PBT_KidneyHiSeq_WorkingData_C.Maher/DE-analysis/2_alignedData/trinity-sep/G1_KI_trinity/trinity

grep -c ">" trinity.fasta
46538

cd ../../G2_KI_trinity/trinity/
grep -c ">" trinity.fasta
43283

cd ../../G3_KI_trinity/trinity/
grep -c ">" trinity.fasta
35991

cd ../../G4_KI_trinity/trinity/
grep -c ">" trinity.fasta
39239

cd ../../G5_KI_trinity/trinity/
grep -c ">" trinity.fasta
37379

cd ../../G6K_CDNJBANXX_trinity/trinity/
grep -c ">" trinity.fasta
214799

cd ../../G7K_CDNJBANXX_trinity/trinity/
grep -c ">" trinity.fasta
266539

cd ../../G8_KI_trinity/trinity/
grep -c ">" trinity.fasta
43050

cd ../../../trinity-all/trinity/
grep -c ">" trinity.fasta
373287









