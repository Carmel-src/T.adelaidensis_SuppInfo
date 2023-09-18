Carmel@CarmelASUS MINGW64 ~
$ cd "E:/@DATA/@@PBT_KidneyIsoSeq_Working Data_C.Maher/IsoSeq-analysis/data/Seqtk/reference COPIES 2020-11"

2020-11-17
-----

cd "C:\Users\Carmel\Desktop\reference COPIES 2020-11"

for F in *.*
	do
		wc -l "$F"
	done


for F in *.*
	do
		grep -C ">" "$F"
	done

-----

19814 clstrTranscriptID.fa					/2 = 9907
9907 clstrTranscriptID.lst

hq.fasta.no5merge.collapsed.rep.fa
pygmy.ANGEL.cds
pygmy.ANGEL.pep
pygmy.ANGEL.utr

-----
# Sort and filter clstrTranscriptID.lst by nuique values to make sure there are no dupicates
sort clstrTranscriptID.lst | uniq > clstrUniq.txt
wc -l clstrUniq.txt
9907 clstrUniq.txt

#List duplicates if there are any found (there are not)
sort clstrTranscriptID.lst | uniq -d > clstrUniqD.txt
wc -l clstrUniqD.txt
0 clstrUniqD.txt


**********


Carmel@CarmelASUS MINGW64 /e/@DATA/@@PBT_KidneyIsoSeq_Working Data_C.Maher/IsoSeq-analysis/data/Seqtk/reference COPIES 2020-11


# ---> 1. Isoseq HQ output file: hq.fasta.no5merge.collapsed.rep.fa


# header format
head -n 2 hq.fasta.no5merge.collapsed.rep.fa
>PB.1.1|000643|path0:1-1879(+)|transcript/14742 transcript/14742 full_length_coverage=3;length=1887;num_subreads=52
ATAGAGTACAAGTACTAGAGTTTTAGTGTGGTCGAAGAGAAATCCAAGATCCACCATGGG

# no of sequences
grep -c ">" hq.fasta.no5merge.collapsed.rep.fa
15729


# ---> 2. Angel CDS output from Isoseq HQ file: pygmy.ANGEL.cds

# header format
head -n 2 pygmy.ANGEL.cds
>PB.2.1|002537|path0:1-1624(+)|transcript/18304|m.1 type:likely-NA len:135 strand:+ pos:282-686
ATGTCCTACTGCAGGCAAGAAGGAAAGGACAGAATTATATTTGTGACAAAAGAGGACCATGAAACTCCAAGTAGTGCTGAGTTAGTAGCGGATGACCCCAATGATCCTTATGAAGATCAAGGTTTGATATTGCCTAATGGAGATATCAATTGGAATTGTCCTTGTCTCGGTGGAATGGCTAGTGGCCCTTGTGGGGAACAATTCAAATCGGCCTTTTCCTGTTTCCATTATAGCAAGGAAGAAATAAAGGGATCAGATTGTGTGGACCAGTTTCGTGCAATGCAAGAGTGCATGCAGAAGTACCCAGATCTATACCCCCAAGAGGAAGATGATGAAGAAGAAAAGCAGAGCAAAGAATTAGAAACAGCTCCTTCAGTAGCCAAAGAGGAAGAGGGATCTAGCTAA

# Note Sequence is on a single line

# no of sequences
grep -c ">" pygmy.ANGEL.cds
13882


# ---> Angel CDS file clusterd with cdHit: 1598936109.fas.1.clstr 

# header format
cd cd-hit\ 99
head 1598936109.fas.1.clstr
>Cluster 0
0       2513aa, >PB.715.4|1b393a|path6:6-8031(+)|transcript/8|m.1216... *
1       2250aa, >PB.715.5|1b393a|path6:278-8030(+)|transcript/13|m.1217... at 99.56%
2       2369aa, >PB.715.6|1b393a|path6:524-8029(+)|transcript/20|m.1218... at 99.58%
3       2148aa, >PB.715.7|1b393a|path6:1187-8030(+)|transcript/49|m.1219... at 99.53%
4       1913aa, >PB.715.9|1b393a|path6:1893-8028(+)|transcript/118|m.1221... at 99.79%
5       1688aa, >PB.715.11|1b393a|path6:2596-8054(+)|transcript/279|m.1223... at 99.59%
6       1452aa, >PB.715.14|1b393a|path6:3304-8029(+)|transcript/591|m.1226... at 100.00%
7       1333aa, >PB.715.15|1b393a|path6:3629-8031(+)|transcript/935|m.1227... at 99.77%
8       1314aa, >PB.715.16|1b393a|path6:3916-8078(+)|transcript/962|m.1228... at 100.00%

# File contains truncated headers from the Angel CDS file, clustered acoording to similarity. Truncation appears to be at the first whitespace

# no of sequences
grep -c ', >' 1598936109.fas.1.clstr     #search pattern changed to avoid the ">" in front of each cluster label eg >Cluster 0
13882
# Sequence number matches Angel CDS file

# no of clusters
grep -c '>C' 1598936109.fas.1.clstr
9907

grep -c '>C' 1598936109.fas.1.clstr.sorted
9907
#both match

tail 1598936109.fas.1.clstr
>Cluster 9902
0       99aa, >PB.8346.1|transcript/22300:1-1143(+)|transcript/22300|m.14112... *
>Cluster 9903
0       99aa, >PB.8372.1|transcript/22508:1-1103(+)|transcript/22508|m.14137... *
>Cluster 9904
0       99aa, >PB.8489.1|transcript/23315:1-855(+)|transcript/23315|m.14244... *
>Cluster 9905
0       99aa, >PB.8706.1|transcript/25054:1-372(+)|transcript/25054|m.14469... *
>Cluster 9906
0       99aa, >PB.9028.1|transcript/4242:1-2955(+)|transcript/4242|m.14776... *

tail 1598936109.fas.1.clstr.sorted
>Cluster 9902
0       99aa, >PB.8346.1|transcript/22300:1-1143(+)|transcript/22300|m.14112... *
>Cluster 9903
0       99aa, >PB.8372.1|transcript/22508:1-1103(+)|transcript/22508|m.14137... *
>Cluster 9904
0       99aa, >PB.8489.1|transcript/23315:1-855(+)|transcript/23315|m.14244... *
>Cluster 9905
0       99aa, >PB.8706.1|transcript/25054:1-372(+)|transcript/25054|m.14469... *
>Cluster 9906
0       99aa, >PB.9028.1|transcript/4242:1-2955(+)|transcript/4242|m.14776... *

# Cluster numbers also match both files:Cluster 9906 + Cluster 0 = 9907

# How has the header changed from Isoseq HQ after processing with Angel
cd ../
grep '>PB.715.4|1b393a|path6:6-8031(+)|transcript/8' hq.fasta.no5merge.collapsed.rep.fa
>PB.715.4|1b393a|path6:6-8031(+)|transcript/8 transcript/8 full_length_coverage=43;length=7837;num_subreads=60

# the "|m.xxxxx" suffix is added after the first whitespace


# ---> 3. The longest seed sequence from the cdHit output: clstrUniq.txt

# header format
head clstrUniq.txt
PB.10.1|004815|path1:5-1713(+)|transcript/17534|m.3
PB.100.1|049515|path5:1-1563(+)|transcript/19243|m.132
PB.1000.1|26d4b2|path0:1-2889(+)|transcript/4521|m.1797
PB.1000.2|26d4b2|path0:1287-2889(+)|transcript/18232|m.1798
PB.1002.1|26ea71|path2:1-2673(+)|transcript/6271|m.1799
PB.1002.2|26ea71|path2:11-2691(+)|transcript/8089|m.1800
PB.1003.1|26fc93|path2:1-3987(+)|transcript/1197|m.1801
PB.1005.1|2718e2|path0:1-5294(+)|transcript/258|m.1804
PB.1006.1|272202|path0:1-2713(+)|transcript/5814|m.1805
PB.1007.1|272202|path1:1-3370(+)|transcript/2554|m.1806

# Prefix ">" and suffix "... *" have been removed in R Studio

# no of strings
grep -c 'PB.' clstrUniq.txt
9907

#matches number of clusters above


# strip the "|m.xxxx" from the strings in clstrUniq.txt 
sed 's/|m.*//' clstrUniq.txt > new_cluster.txt

# check
head new_cluster.txt
PB.10.1|004815|path1:5-1713(+)|transcript/17534
PB.100.1|049515|path5:1-1563(+)|transcript/19243
PB.1000.1|26d4b2|path0:1-2889(+)|transcript/4521
PB.1000.2|26d4b2|path0:1287-2889(+)|transcript/18232
PB.1002.1|26ea71|path2:1-2673(+)|transcript/6271
PB.1002.2|26ea71|path2:11-2691(+)|transcript/8089
PB.1003.1|26fc93|path2:1-3987(+)|transcript/1197
PB.1005.1|2718e2|path0:1-5294(+)|transcript/258
PB.1006.1|272202|path0:1-2713(+)|transcript/5814
PB.1007.1|272202|path1:1-3370(+)|transcript/2554

grep -c 'PB.' new_cluster.txt
9907

# How many of thes strings are in the original ONE LINE isoseq file
grep -c -Ff new_cluster.txt hq.fasta.no5merge.collapsed.rep.1L.fa 
9813

# 94 strings are missing


# ---> 4. find the missing strings

# convert the isoform hq file into the same format as new_cluster.txt
grep '>' hq.fasta.no5merge.collapsed.rep.1L.fa| sed s'/ tra.*//' | sed s'/>//' > hq_headers 
head hq_headers 
PB.1.1|000643|path0:1-1879(+)|transcript/14742
PB.2.1|002537|path0:1-1624(+)|transcript/18304
PB.3.2|00361c|path0:43-2240(+)|transcript/10564
PB.4.1|004079|path1:6-745(+)|transcript/23781
PB.4.2|004079|path1:6-612(+)|transcript/24311
PB.5.1|004079|path2:209-642(+)|transcript/24981
PB.6.1|004079|path3:273-616(+)|transcript/25193
PB.7.1|004079|path5:6-664(+)|transcript/24242
PB.8.1|004079|path6:1-607(+)|transcript/24477
PB.9.1|004079|path7:1-725(+)|transcript/23873

cat hq_headers |wc -l
15729

# just checking I still get the same number of matches
grep -Ff new_cluster.txt hq_headers > matching_strings
cat matching_strings |wc -l
9813

#check for duplicates
cat hq_headers | sort |uniq -cd
0

cat new_cluster.txt | sort |uniq -cd
      2 PB.1206.1|2e8181|path5:1-4617(+)|transcript/551
      2 PB.1261.1|30ba06|path3:633-1922(+)|transcript/21324
      2 PB.1570.2|3bd5a8|path2:20-1934(+)|transcript/13977
      2 PB.1719.1|404d72|path0:1-3429(+)|transcript/2407
      2 PB.1764.3|429969|path8:44-2062(+)|transcript/12739
      2 PB.1813.1|446aed|path0:1-2871(+)|transcript/4320
      2 PB.1865.1|46bd4e|path0:1-2536(+)|transcript/6726
      2 PB.1869.1|46d663|path2:1-1320(+)|transcript/21375
      2 PB.220.1|08884f|path0:1-4120(+)|transcript/991
      2 PB.2360.1|5843a6|path32:1-2062(+)|transcript/11326
      2 PB.2541.3|5e93d1|path0:5-2495(+)|transcript/7381
      2 PB.2560.1|5f09e2|path0:1-1196(+)|transcript/21957
      2 PB.2618.1|612418|path1:1-2013(+)|transcript/12200
      2 PB.2713.3|65abaa|path0:37-2495(+)|transcript/7693
      2 PB.2892.1|6b359a|path12:1-4597(+)|transcript/526
      2 PB.3037.2|702e9b|path5:2-3018(+)|transcript/3923
      2 PB.3243.1|784f79|path0:1-2683(+)|transcript/5903
      2 PB.3259.1|792211|path3:1-5233(+)|transcript/268
      2 PB.3286.3|7a2645|path3:103-3087(+)|transcript/3780
      2 PB.3405.1|7f2f96|path0:5-2081(+)|transcript/11903
      3 PB.3430.1|8055ae|path8:1-3589(+)|transcript/2043
      2 PB.3461.1|81f072|path5:1-2492(+)|transcript/7285
      2 PB.351.4|0cccb8|path5:27-4293(+)|transcript/775
      2 PB.3562.1|864d1b|path0:1-3492(+)|transcript/2174
      2 PB.3615.1|883e51|path3:1-4044(+)|transcript/1347
      2 PB.3836.3|8ff0e9|path0:4-3972(+)|transcript/946
      2 PB.3934.1|93b547|path0:1-3353(+)|transcript/2676
      2 PB.4042.2|97fba8|path2:34-2096(+)|transcript/13603
      2 PB.4219.2|9d8e09|path3:9-3566(+)|transcript/2045
      2 PB.4442.1|a6620e|path0:1-2906(+)|transcript/4555
      2 PB.4458.1|a6e2a5|path0:1-3673(+)|transcript/1739
      2 PB.4822.1|b4950f|path1:1-3315(+)|transcript/2177
      2 PB.4822.2|b4950f|path1:1240-3315(+)|transcript/10982
      2 PB.4919.1|b7da5b|path0:1-3327(+)|transcript/2652
      2 PB.4923.1|b7fed9|path0:1-2394(+)|transcript/8643
      2 PB.5080.1|bde542|path3:1-3060(+)|transcript/3777
      2 PB.5096.2|be8c26|path4:70-2754(+)|transcript/5782
      2 PB.5149.1|c15864|path0:1-1628(+)|transcript/17934
      2 PB.5308.1|c8b5b5|path2:35-3805(+)|transcript/1487
      2 PB.5395.1|cbd93f|path0:1-1699(+)|transcript/17107
      2 PB.5483.2|cf1828|path1:23-2344(+)|transcript/9126
      2 PB.5619.1|d4d2ea|path0:1-2677(+)|transcript/6673
      2 PB.5663.1|d5e71b|path0:1-5101(+)|transcript/303
      2 PB.5939.1|e23f0e|path0:1-2960(+)|transcript/3907
      2 PB.5939.2|e23f0e|path0:5-4542(+)|transcript/622
      2 PB.613.1|16d05c|path0:1-1558(+)|transcript/18673
      2 PB.6133.1|ea032e|path9:1-2052(+)|transcript/12270
      2 PB.6276.1|ee992b|path0:1-1487(+)|transcript/19762
      2 PB.6283.1|eee48e|path1:1-3454(+)|transcript/2331
      2 PB.6286.1|ef0b10|path12:1-2744(+)|transcript/5222
      2 PB.6303.1|efc514|path0:1-3222(+)|transcript/2889
      2 PB.6378.3|f31894|path9:80-2754(+)|transcript/5975
      2 PB.6385.1|f32c31|path0:1-1500(+)|transcript/12840
      2 PB.6426.1|f4ded6|path0:1-1796(+)|transcript/16785
      2 PB.6576.2|facceb|path2:1-3314(+)|transcript/3157
      2 PB.659.1|19437d|path1:1-3483(+)|transcript/2210
      2 PB.6763.1|transcript/10186:1-2230(+)|transcript/10186
      2 PB.6877.1|transcript/1099:1-4113(+)|transcript/1099
      2 PB.6966.1|transcript/11636:1-2103(+)|transcript/11636
      2 PB.713.1|1b043a|path0:1-965(+)|transcript/23154
      2 PB.7227.1|transcript/138:1-5869(+)|transcript/138
      2 PB.7287.1|transcript/14227:1-1898(+)|transcript/14227
      3 PB.7342.1|transcript/14710:1-1879(+)|transcript/14710
      2 PB.7349.1|transcript/14730:1-1903(+)|transcript/14730
      2 PB.7353.1|transcript/14772:1-1860(+)|transcript/14772
      2 PB.74.1|02ebe4|path0:1-2606(+)|transcript/6537
      2 PB.7501.1|transcript/1596:1-3759(+)|transcript/1596
      2 PB.7528.1|transcript/16151:1-1797(+)|transcript/16151
      2 PB.7678.1|transcript/1728:1-3672(+)|transcript/1728
      2 PB.7707.1|transcript/17496:1-1683(+)|transcript/17496
      2 PB.7868.1|transcript/188:1-5558(+)|transcript/188
      2 PB.7973.1|transcript/19628:1-1476(+)|transcript/19628
      2 PB.8082.1|transcript/20433:1-1412(+)|transcript/20433
      2 PB.816.1|1e334e|path4:1-2040(+)|transcript/12316
      2 PB.8325.1|transcript/22158:1-1191(+)|transcript/22158
      2 PB.834.2|1ec273|path2:1291-2750(+)|transcript/19945
      2 PB.8389.1|transcript/22625:1-1060(+)|transcript/22625
      2 PB.8709.1|transcript/252:1-5299(+)|transcript/252
      2 PB.8883.1|transcript/3393:1-3152(+)|transcript/3393
      2 PB.8886.1|transcript/3402:1-3146(+)|transcript/3402
      2 PB.8920.1|transcript/3588:1-3076(+)|transcript/3588
      2 PB.9068.1|transcript/4442:1-2918(+)|transcript/4442
      2 PB.9162.1|transcript/4983:1-2815(+)|transcript/4983
      2 PB.9209.1|transcript/5238:1-2771(+)|transcript/5238
      2 PB.924.1|22d3cf|path3:1-1494(+)|transcript/19025
      2 PB.9277.1|transcript/5650:1-2707(+)|transcript/5650
      2 PB.944.1|241290|path4:2-2629(+)|transcript/6557
      2 PB.9535.1|transcript/7343:1-2508(+)|transcript/7343
      2 PB.9581.1|transcript/7647:1-2463(+)|transcript/7647
      2 PB.9683.1|transcript/8355:1-2402(+)|transcript/8355
      2 PB.9766.1|transcript/9024:1-2324(+)|transcript/9024
      2 PB.9784.1|transcript/915:1-4197(+)|transcript/915

# OK so during the clustering it appears that several of the sequences have been made seed sequences (denoted with a *) for more than one cluster.

# The counts above equate to how many strings
cat new_cluster.txt | sort |uniq -cd | cut -d' ' -f7 | awk '{sum+=$1-1} END{printf("%d\n",sum)}'
94

# Which is the discrepancy noted in section 3 above and using the non manipulated multi line fasta file. 


# Final pulling out of transcripts on a single line 
# note: NOT using seqtk anymore:

grep -Ff new_cluster.txt -A 1 hq.fasta.no5merge.collapsed.rep.1L.fa > reference_transcripts.1L.fasta

grep -c ">" reference_transcripts.1L.fasta
9813


# also remove the lines that grep outputs as '--' as this upsets bbduk

grep -Ff new_cluster.txt -A 1 hq.fasta.no5merge.collapsed.rep.1L.fa | grep -v ^-- > reference_transcripts.1Lv.fasta

grep -c ">" reference_transcripts.1Lv.fasta
9813


# -----next step------

# due to product/DeepThought availability run bbduk on all 3 files:

reference_transcripts.1Lv.fasta
reference_transcripts_seqtk.1L.fasta
reference_transcripts_seqtk.fasta

# remove-poly-a_bbduk.bash ran on all 3 using 2 program methods 20/12/2020

# however reference_transcripts.1L.fasta (or reference_transcripts.1Lv.fasta) is the master file assuming all else is equal
# This trimmed file will then be used as the reference for Kallisto on NECTAR

# Trim polya command in bbduk leaves some tails intact due to single bases and does not facilitate hamming distance command at the same time (?), however specifying Ktrim=r literal or "twice_trimmed" (poly-a and THEN ktrim) sequences are very short.
# tried with different K values. 

# Terry Bertozzi provided perl script as below to manually trim.


# -----next step------

# Terry Bertozzi
#command line to remove poly A tails

#for use on /mnt/IsoSeq-analysis/data/Seqtk/reference_transcripts.1Lv.fasta

cd /mnt/IsoSeq-analysis/data/Seqtk

perl -pe 's/(A{5,}.{0,2})+$//gm' reference_transcripts.1L.fasta > reference_transcripts.1L.clean.fasta

perl -pe 's/(A{5,}.{0,2})+$//gm' reference_transcripts.1Lv.fasta > reference_transcripts.1Lv.clean.fasta

grep -c ">" reference_transcripts.1Lv.clean.fasta
9813

head reference_transcripts.1Lv.clean.fasta
#also looks good.


reference_transcripts.1Lv.clean.fasta
#IS THE REFERENCE


This command now works (worked on reference_transcripts.1L.fasta) re-run on the 1Lv file that does not have the -- grep line placeholders 

