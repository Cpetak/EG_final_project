---
title: "Mapping Workflow"
author: "Baxter Worthing"
date: "4/15/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

The exome data are already mapped to the reduced reference, but the 77 sets of RNA-seq reads are not.

here is my code for mapping them 

```{bash, eval =F}

# navigate to our project dir
cd /data/project_data/GroupProjects/UTR

# grab the reduced reference
cp /data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa referenceSeqs/

# make bowtie index with the reference 

bowtie2-build referenceSeqs/Pabies1.0-genome_reduced.fa Pabies1.0-genome_reduced

# grab annotation 

cp /data/project_data/Annotation/Pabies01-gene.gff3.gz annotations/

# check names of annotation vs index

bowtie-inspect --names Pabies1.0-genome_reduced

# test tophat command with just one fastq file 
# using --transcriptome-index so I can drop -G in the future
# note I already did ASC_06_C_0_GAGTCC_R1 as well 

tophat2 -G /data/project_data/GroupProjects/UTR/annotations/Pabies01-gene.gff3 --transcriptome-index /data/project_data/GroupProjects/UTR/annotations/reducedRefTopHat --num-threads 3 -o /data/project_data/GroupProjects/UTR/thRNAseqBams/ASC_06_C_10_ACGTCT_R1 Pabies1.0-genome_reduced /data/project_data/RS_RNASeq/fastq/cleanreads/ASC_06_C_10_ACGTCT_R1.cl.fq 


```

This command works, it runs in about 20 min and gave mapping rate of ~60% for the two sets of reads I tried it with 

now loop over all the other sets of RNAseq reads with a bash script

```{bash, eval =F}
cd referenceSeqs/

bash ../mapBatchTH.sh

```
