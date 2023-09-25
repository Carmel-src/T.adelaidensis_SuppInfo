# Transcript Annotation and Seasonal Gene Expression in Kidney Tissue of the Australian Pygmy Bluetongue Lizard, _Tiliqua adelaidensis_ 
## Supplementary Information - Included in Thesis Appendix #2

This repository contains the analysis process used for transcript assembly, annotation, and gene expression in kidney tissue of the Australian pygmy bluetongue skink _Tiliqua adelaidensis_.  

### Methods are best followed through the [full workflow document](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/blob/main/Full-Workflow-2023.pdf) rather than browsing uploaded files
This PDF file of the full workflow has bookmarks and an interactive contents page for navigation of analysis sections. This is also provided as the raw [rmd](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/blob/main/Full-Workflow-2023.Rmd) file and a [html](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/blob/main/Full-Workflow-2023.html) output.

This document outlines the full pipeline and links to script files in the context of each step. Most of the data analysis for Annotation & Gene Expression was conducted in R and is presented in the full workflow document, therefore there are fewer additional files in those folders here.


***


The Folders in this repository contain supplementary information including the full script files organised per chapter, which are linked to in the relevant section of the above workflow:

#### [Assembly](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/tree/main/Assembly) - De-novo Assembly of short-read ranscripts

This section of the document outlines _de-novo_ assembly of long-read Iso-seq data sequenced from poly-a selected RNA extracted from a single Kidney. These sequences will be used to annotate genes found in kidney tissue of _T. adelaidensis_ and create a set of reference transcripts for gene expression analysis conducted in later chapters.  
This section also outlines _de-novo_ assembly of short-reads sequenced on an Illumina HiSeq and compares these assemblies.  


#### [Annotation](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/tree/main/Annotation) - Annotation of Transcripts

This section of the document outlines annotation of transcripts (and isoforms) of the full length long read transcripts as assembled in Chapter 4.  


#### [Gene Expression](https://github.com/Carmel-src/T.adelaidensis_SuppInfo/tree/main/Gene%20Expression) - Gene Expression Analysis between Seasons

This section of the document outlines gene expression analysis conducted in Chapter 5, of the same short-read data presented in Chapter 3 supplementary information, however instead of using _de-novo_ assemblies, cleaned reads were directly aligned to the long read reference produced in the previous sections. Much of this analysis was then conducted in R Studio so this chapter especially, is best followed through the Appendix document.
