# Our Lab Nootbook for Ecological Genomics final project

## Author: Baxter Worthing, Csenge Petak, Kerria Burns

### Start Date: 2020-04-15
### End Date: 2020-05-015
### Project Descriptions:   
Final project lab book


# Table of Contents:   
* [Entry 1: 2020-04-15, Wednesday](#id-section1)


------    
<div id='id-section1'/>   


### Entry 1: 2020-04-15, Wednesday.   
These are all the steps we took / are planning to take. Order doesn't always matter.

####Step 0: Building index for reduced reference

bowtie2-build referenceSeqs/Pabies1.0-genome_reduced.fa Pabies1.0-genome_reduced

Step 1: Mapping transcriptomic data to reduced reference with Tophat

Step 2: Using samtools to create vcf

Step 3: Installing SNPeff

Step 4. Analysing phenotypic data in R using ANOVA and ggplot

Step 5. Using bioclim data to find more pops in HD and CW

Step 6. Creating genotype likelihood for all individuals in CW and all individuals in HD

/data/project_data/GroupProjects/UTR/ANGSD_GL

HD bam files = HD_bamlist.list
CW bam files = CW_bamlist.list

Script for all sites (including monomorphic), HD, folded
This was used for GL that would be later used for nucleotide diversity analysis
``````
#!/bin/bash

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"
output="/data/project_data/GroupProjects/UTR/ANGSD_GL/folded_GL_all"

# Estimating GL's and allele frequencies for all sites with ANGSD

ANGSD -b /data/project_data/GroupProjects/UTR/ANGSD_GL/HD_bamlist.list \
-ref ${REF} -anc ${REF} \
-out ${output}/HD_allsites \
-nThreads 2 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
-fold 1
# -SNP_pval 1e-6, this is uncommented in GL prepared for Fst where we want only polymorhic sites

``````

Step 7. Comparing global nucleotide diversities between CW and HD

One of the outputs of the previous script is then used for "realSFS" as below:
``````
realSFS /data/project_data/GroupProjects/UTR/ANGSD_GL/test_results/test01_allsites.saf.idx -maxIter 1000 -tole 1e-6 -P 1 > /data/project_data/GroupProjects/UTR/ANGSD_GL/test_results/test01_allsites.sfs
``````
test01 is replaced with HD or CW accordingly

Then, thetas are calculated as such:
``````
#!/bin/bash

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"
output="/data/project_data/GroupProjects/UTR/ANGSD_GL/test_results"

# Estimating thetas for all sites, using the SFS from above as a prior to estimate the GL's
ANGSD -b /data/project_data/GroupProjects/UTR/ANGSD_GL/test2list.list \
-ref ${REF} -anc ${REF} \
-out ${output}/test03_allsites \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-setMinDepth 3 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doSaf 1 \
-pest ${output}/test03_allsites.sfs \
-doThetas 1 \
-fold 1

thetaStat do_stat ${output}/test03_allsites.thetas.idx
``````

Step 8. Getting global Fst between CW and HD

``````
#!/bin/bash

output="/data/project_data/GroupProjects/UTR/ANGSD_GL/test_results"

realSFS ${output}/test01_allsites.saf.idx ${output}/test02_allsites.saf.idx -P 1 > ${output}/test01_test02_allsites.sfs
realSFS fst index ${output}/test01_allsites.saf.idx ${output}/test02_allsites.saf.idx -sfs ${output}/test01_test02_allsites.sfs -fstout ${output}/test01_test02_allsites -whichFst 1
realSFS fst print ${output}/test01_test02_allsites.fst.idx > ${output}/test01_test02_allsites.fst
``````
replace test01 with HD, test02 with CW


Step 9. Identifying UTRs

Step 10. Getting Fst values for single SNPs

Step 11. Getting per base nuvleotide diversities

Step 12. Check for functional enrichment for high Fst SNPs

Step 13. Testing for selection at the per-SNP level using PCAngsd


------    
