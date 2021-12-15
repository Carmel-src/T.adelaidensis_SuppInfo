
CDHit notes:


http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=Server%20home

cd-hit
CD-HIT clusters proteins that meet a similarity threshold, usually a sequence identity. Each cluster has one representative sequence. The input is a protein dataset in fasta format. It generates a fasta file of representative sequences and a text file of list of clusters.

cd-hit-est
CD-HIT-EST clusters a nucleotide sequences that meet a similarity threshold, usually a sequence identity. The input is a DNA/RNA dataset in fasta format It generates a fasta file of representative sequences and a text file of list of clusters. It can not be used for very long sequences, like full genomes.

h-cd-hit
Multiple CD-HIT runs. Proteins are first clustered at a high identity (like 90%), the non-redundant sequences are further clustered at a low identity (like 60%). A third cluster can be performed at lower identity. Multi-step run is more efficient and more accurate than a single run.


INput files:
Nucleotides:  E:\@DATA\@@PBT_KidneyIsoSeq_Working Data_C.Maher\IsoSeq-analysis\data\ANGEL\pygmy.ANGEL.cds
Proteins:  E:\@DATA\@@PBT_KidneyIsoSeq_Working Data_C.Maher\IsoSeq-analysis\data\ANGEL\pygmy.ANGEL.pep


cd-hit - run on pygmy.ANGEL.pep


separate run cutoff at .99

***
Your job 1598936109 is finished.
Program you ran: cd-hit
You input file is pygmy.ANGEL.pep and we named it as 1598936109.fas.0
Summary information for 1598936109.fas.0 included in 1598936109.fas.0.stat
You required 1 runs for sequence clustering
     1. Fasta file for representative sequences at 99% identity is 1598936109.fas.1
         Summary information for 1598936109.fas.1 included in 1598936109.fas.1.stat
         Corresponding cluster file is1598936109.fas.1.clstr
         Sorted cluster file by size is 1598936109.fas.1.clstr.sorted
Generated shell script is run-1598936109.sh

faa_stat.pl 1598936109.fas.0
faa_stat.pl 1598936109.fas.1
/data5/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1598936109.fas.1.clstr > 1598936109.fas.1.clstr.sorted
/data5/data/NGS-ann-project/apps/cd-hit/clstr_list.pl 1598936109.fas.1.clstr 1598936109.clstr.dump
gnuplot1.pl < 1598936109.fas.1.clstr > 1598936109.fas.1.clstr.1; gnuplot2.pl 1598936109.fas.1.clstr.1 1598936109.fas.1.clstr.1.png
/data5/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1598936109.clstr.dump 1598936109.clstr_no.dump
/data5/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1598936109.clstr.dump 1598936109.clstr_len.dump len
/data5/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1598936109.clstr.dump 1598936109.clstr_des.dump des
***



You need to unzip within git bash - further cluster manipulation within R

cd /e/@DATA/@@PBT_KidneyIsoSeq_Working Data_C.Maher/IsoSeq-analysis/data/cdhit/pep

tar -xvf <file>
 
 