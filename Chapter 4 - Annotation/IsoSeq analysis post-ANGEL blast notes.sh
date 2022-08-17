
#Wherever you run this script is where the .out will appear.
/mnt/Prog/blast/blastdb/AcarProt
/mnt/Prog/blast/blastdb/sprot

BLAST 2.9.0+
Uniprot sprot & trembl databases downloaded: 28/5/19
Anolis .pep fasta files downloaded 12/11/19

These lines are run from within a screen on the NECTAR machine. Full paths are included to avoid errors in working directory placement.


******

##this group is named as: pygmy.ANGEL_blastx_<db>.out


#BLAST the pygmy.ANGEL.cds against the protein database
#create database from within db folder
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/makeblastdb -in uniprot_sprot.fasta -parse_seqids -blastdb_version 5 -title "sprot" -dbtype prot -out sprot

#search
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx –db /mnt/Prog/blast/blastdb/sprot –query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds –out pygmy.ANGEL_blastx_uniprot.out -outfmt 6 -evalue 0.00001



#BLAST the pygmy.ANGEL.cds against the anolis proteome
#create database from within db folder

#AnoCar2.0.pep.all
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/makeblastdb -in Anolis_carolinensis.AnoCar2.0.pep.all.fa -parse_seqids -blastdb_version 5 -title "sprot" -dbtype prot -out AnoCar2.0.pep.all

#search
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all.out -outfmt 6 -evalue 0.00001

#AnoCar2.0.pep.abinitio
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/makeblastdb -in Anolis_carolinensis.AnoCar2.0.pep.abinitio.fa -parse_seqids -blastdb_version 5 -title "sprot" -dbtype prot -out AnoCar2.0.pep.abinitio

#search
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.abinitio/AnoCar2.0.pep.abinitio -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.abinitio.out -outfmt 6 -evalue 0.00001


******

##this group is named as: pygmy.ANGEL_blastx_<db>-maxtarg5-maxhsp5.out

# uniprot sprot
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/sprot/sprot -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_sprot-maxtarg5-maxhsp5.out -outfmt 6 -max_target_seqs 5 -max_hsps 5 -evalue 0.00001 -num_threads 3

#AnoCar2.0.pep.all
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all-maxtarg5-maxhsp5.out -outfmt 6 -max_target_seqs 5 -max_hsps 5 -evalue 0.00001 -num_threads 3

#AnoCar2.0.pep.abinitio
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.abinitio/AnoCar2.0.pep.abinitio -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.abinitio-maxtarg5-maxhsp5.out -outfmt 6 -max_target_seqs 5 -max_hsps 5 -evalue 0.00001 -num_threads 1


******

##this group is named as: pygmy.ANGEL_blastx_<db>-eval1e-10.out

# uniprot sprot
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/sprot/sprot -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_sprot-eval1e-10.out -outfmt 6 -evalue 1e-10 -num_threads 3

#AnoCar2.0.pep.all
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all-eval1e-10.out -outfmt 6 -evalue 1e-10 -num_threads 3

#AnoCar2.0.pep.abinitio
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.abinitio/AnoCar2.0.pep.abinitio -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.abinitio-eval1e-10.out -outfmt 6 -evalue 1e-10 -num_threads 1


******

##this group is named as: pygmy.ANGEL_blastx_<db>-maxtarg1-maxhsp5.out
Warning: [blastx] Examining 5 or more matches is recommended

# uniprot sprot
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/sprot/sprot -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_sprot-maxtarg1-maxhsp5.out -outfmt 6 -max_target_seqs 1 -max_hsps 5 -evalue 0.00001 -num_threads 1

#AnoCar2.0.pep.all
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all-maxtarg1-maxhsp5.out -outfmt 6 -max_target_seqs 1 -max_hsps 5 -evalue 0.00001 -num_threads 3

#AnoCar2.0.pep.abinitio
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.abinitio/AnoCar2.0.pep.abinitio -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.abinitio-maxtarg1-maxhsp5.out -outfmt 6 -max_target_seqs 1 -max_hsps 5 -evalue 0.00001 -num_threads 3


******

#Format 5 for BLAST2GO     ----     unisprot & anolis_all
# -max_target_seqs 5 -max_hsps 5 -evalue 0.00001

/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all-BLAST2GO -outfmt 5 -max_target_seqs 5 -max_hsps 5 -evalue 0.00001 -num_threads 4

/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/sprot/sprot -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_sprot-BLAST2GO -outfmt 5 -max_target_seqs 5 -max_hsps 5 -evalue 0.00001 -num_threads 4


******

##this group is named as: pygmy.ANGEL_blastx_<db>-maxtarg1-maxhsp1.out
Warning: [blastx] Examining 5 or more matches is recommended

# uniprot sprot
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/sprot/sprot -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_sprot-maxtarg1-maxhsp1.out -outfmt 6 -max_target_seqs 1 -max_hsps 1 -evalue 0.00001 -num_threads 3

#AnoCar2.0.pep.all
/mnt/Prog/blast/ncbi-blast-2.9.0+/bin/blastx -db /mnt/Prog/blast/blastdb/AcarProt/AnoCar2.0.pep.all/AnoCar2.0.pep.all -query /mnt/IsoSeq-analysis/data/ANGEL/pygmy.ANGEL.cds -out pygmy.ANGEL_blastx_AnoCar.pep.all-maxtarg1-maxhsp1.out -outfmt 6 -max_target_seqs 1 -max_hsps 1 -evalue 0.00001 -num_threads 3
