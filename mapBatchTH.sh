#!/bin/bash

# map all the RNA seq reads to the reduced reference genome using Tophat2 and default params
# needs to be run in dir where bowtie index is will send output to /users/b/w/bwworthi/groupProj/mapping/cds2kbBams

# gab all the paths to all the fastq files
# need to ignore that one blank file and the 2 I already did

for f in `find /data/project_data/RS_RNASeq/fastq/cleanreads/ -type f| grep -E -v "BLANK|ASC_06_C_0_GAGTCC_R1|ASC_06_C_10_ACGTCT_R1"`

do


# trim down name for output naming
outName=`basename $f | cut -f1 -d'.'`


# outname looks like XBM_07_H_10_GGAGGT_R1

#mkdir /data/project_data/GroupProjects/UTR/thRNAseqBams/${outName}


# run tophat2 to map to the reduced P. abies reference with reads designated as $f (at the end) and send to its own output dir
# uses pre-made transcriptome index to save time  
tophat2  --transcriptome-index /data/project_data/GroupProjects/UTR/annotations/reducedRefTopHat --num-threads 3 \
-o /data/project_data/GroupProjects/UTR/thRNAseqBams/${outName} Pabies1.0-genome_reduced $f

done
