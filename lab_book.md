# Our Lab Nootbook for Ecological Genomics final project

## Author: Baxter Worthing, Csenge Petak, Kerria Burns

### Start Date: 2020-04-15
### End Date: 2020-05-015
### Project Descriptions:   
Final project lab book


------    
  


### Methods.   
These are all the steps we took / are planning to take. Order doesn't always matter.

#### Step 0: Building index for reduced reference

bowtie2-build referenceSeqs/Pabies1.0-genome_reduced.fa Pabies1.0-genome_reduced

#### Step 1: Mapping transcriptomic data to reduced reference with Tophat
```
# test tophat command with just one fastq file 
# using --transcriptome-index so I can drop -G in the future
# note I already did ASC_06_C_0_GAGTCC_R1 as well 

tophat2 -G /data/project_data/GroupProjects/UTR/annotations/Pabies01-gene.gff3 --transcriptome-index /data/project_data/GroupProjects/UTR/annotations/reducedRefTopHat --num-threads 3 -o /data/project_data/GroupProjects/UTR/thRNAseqBams/ASC_06_C_10_ACGTCT_R1 Pabies1.0-genome_reduced /data/project_data/RS_RNASeq/fastq/cleanreads/ASC_06_C_10_ACGTCT_R1.cl.fq 

# map all 76 fastq files with  2 scripts (saved in this repo)

cd referenceSeqs/

bash ../mapBatchTH.sh
bash ../indexBamsTH.sh

# merge BAM files

cd /data/project_data/GroupProjects/UTR

samtools merge -b RNAbamPaths merged_RNA.bam
```
#### Step 2: Using samtools to create vcf

```
cd /data/project_data/GroupProjects/UTR

# need list of just the BAMs we want for the full analysis 

tail -n +2 HD_CW_extended_list.txt | cut -f1 > ourExomeIDs

find /data/project_data/RS_ExomeSeq/mapping/all/ | grep -f ourExomeIDs | grep -v '.bai$' > ourExomePaths 

# make vcf for all these BAMs

bcftools mpileup --max-depth 1000 -Ou --threads 3 -f /data/project_data/GroupProjects/UTR/referenceSeqs/Pabies1.0-genome_reduced.fa -b /data/project_data/GroupProjects/UTR/ourExomePaths | bcftools call -v -m -Ob -o /data/project_data/GroupProjects/UTR/VCFs/CW_HD_SNPs.bcf

# filter for call quality and convert to vcf 
bcftools view -i '%QUAL>=20' -O z -o VCFs/filtered_CW_HD_SNPs.vcf.gz VCFs/CW_HD_SNPs.bcf

# give it col names
cut -f7 -d'/' ourExomePaths | cut -f1 -d'.'  > HDCWsampleNames

bcftools reheader -s HDCWsampleNames -o VCFs/shortName_filtered_CW_HD_SNPs.vcf.gz VCFs/filtered_CW_HD_SNPs.vcf.gz
```

#### Step 3: Installing SNPeff (we are no longer doing this)

#### Step 4. Analysing phenotypic data in R using ANOVA and ggplot

R script in the phenotypic_data folder in this repo

#### Step 5. Using bioclim data to find more pops in HD and CW

txt file in repo, "RS_Exome_bioclim.txt"

``````
#reading in bioclim data
bioclim <- read.table(file="RS_Exome_bioclim.txt", header=TRUE)
boxplot(bio1_13~Pop, data=bioclim)
boxplot(bio12_13~Pop, data=bioclim, col="blue")


#HotDry (5 fams): ASC_06, BRU_05, ESC_01, XBM_07, NOR_02
#CoolWet (5 fams): CAM_02, JAY_02, KAN_04, LOL_02, MMF_13

library(ggplot2)

ggplot(bioclim, aes(x=bio12_13, y=bio1_13)) + 
  geom_point(aes(color=Pop), show.legend = FALSE) +
  geom_text(label=bioclim$Pop)
``````

#### Step 6. Creating genotype likelihood for all individuals in CW and all individuals in HD

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

getting only polymorphic sites from mafs output

``````
zcat HD_allsites.mafs.gz > HD_allsites.mafs
awk '{if($7<0.000001)print}' <HD_allsites.mafs > test_awk_HD.txt
cat test_awk_HD.txt test_awk_CW.txt > both_poly.txt
sort -k1 -k2 both_poly.txt
uniq -d both_poly.txt #returns only duplicates
wc -l #gives us number of duplicates -> the ones that were polymorphic in both
``````

GL without filtering to get higher resultion per-site nucleotide diversity later on.
``````
#!/bin/bash

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"
output="/data/project_data/GroupProjects/UTR/ANGSD_correct/no_fill"

# Estimating GL's and allele frequencies for all sites with ANGSD

ANGSD -b /data/project_data/GroupProjects/UTR/ANGSD_correct/CW_bamlist.list \
-ref ${REF} -anc ${REF} \
-out ${output}/CW_allsites \
-nThreads 2 \
-GL 1 \
-doSaf 1 \
-fold 1 \
-sites /data/project_data/GroupProjects/UTR/regions_totestG.csv #this allows us to run ANGSD only on specific regions to make it faster. you first have to create a file with the desired locations, format specified here http://www.popgen.dk/angsd/index.php/Sites, then run "angsd sites index your.file", and then refer to this file like I did above when running angsd.
``````

#### Step 7. Comparing global nucleotide diversities between CW and HD

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
R script to visualise SFS

``````
SFS <- scan("HD_allsites.sfs")
sumSFS <- sum(SFS)
pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS_HD <- SFS[2:length(SFS)]
barplot(plotSFS_HD, names.arg = seq(1,length(SFS)-1), xlab="Derived allele frequency", ylab="Number of sites", col=rgb(1, 0, 0, 0.5))

SFS <- scan("CW_allsites.sfs")
sumSFS <- sum(SFS)
pctPoly = 100*(1-(SFS[1]/sumSFS))
plotSFS_CW <- SFS[2:length(SFS)]
barplot(plotSFS_CW, col=rgb(0, 0, 1, .5),add=TRUE)

legend('topright', bty = 'n', title = 'Source climate',
       legend = c('HotDry', 'ColdWet'), fill = c("tomato3", "dodgerblue3"))


div <- read.table("HD_allsites.thetas.idx.pestPG")
colnames(div)=c("window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites
pdf("XWS_diversity_stats2.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
hist(div$tajD,col=rgb(1, 0, 0, 0.5),xlab="Tajima's D",main="")

summary(div)

barplot(plotSFS, main="SFS", xlab= "Derived allele frequency", ylab="Number of sites")
dev.off()

plot(div$tajD)

div <- read.table("CW_allsites.thetas.idx.pestPG")
colnames(div)=c("window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")
div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites
pdf("XWS_diversity_stats2.pdf")
par(mfrow=c(2,2))

hist(div$tWpersite,col="darkorchid",xlab="Theta_W",main="")
hist(div$tPpersite,col="darkorchid", xlab="Theta-Pi",main="")
hist(div$tajD,col=rgb(0, 0, 1, .5),xlab="Tajima's D",main="", add=TRUE)
abline(v = 0, col="red", lwd=3, lty=2)
summary(div)
``````

#### Step 8. Getting global (sliding window?) Fst between CW and HD

``````
#!/bin/bash

output="/data/project_data/GroupProjects/UTR/ANGSD_GL/test_results"

realSFS ${output}/test01_allsites.saf.idx ${output}/test02_allsites.saf.idx -P 1 > ${output}/test01_test02_allsites.sfs
realSFS fst index ${output}/test01_allsites.saf.idx ${output}/test02_allsites.saf.idx -sfs ${output}/test01_test02_allsites.sfs -fstout ${output}/test01_test02_allsites -whichFst 1
realSFS fst print ${output}/test01_test02_allsites.fst.idx > ${output}/test01_test02_allsites.fst
``````
replace test01 with HD, test02 with CW
not sure if this is useful on top of the per site Fst measures

#### Step 9. Identifying UTRs

```
# just pull every CDS out of the gff, turn it into bed and the do bedtools subtract against the reference bed

fgrep 'CDS' annotations/Pabies01-gene.gff3 > annotations/all_CDS.gff3

# remove these from the reference bed

/data/popgen/bedtools2/bin/bedtools subtract -a annotations/reducedFull.bed -b annotations/all_CDS.gff3 > annotations/noCDS_reduced.bed

# use samtools depth with this bed file and the merged bam to see where reads align to hypothetical UTRs

samtools depth -b annotations/noCDS_reduced.bed merged_RNA.bam  > final_RNAbamDepth

# filter out low depth bases 
egrep -v -w '0|1|2|3' final_RNAbamDepth | sed -r 's/[[:space:]]/,/g' > final_RNAbamDepth.csv 

# see R script in this repo for processing of final_RNAbamDepth.csv

```

#### Step 10. Getting Fst values for single SNPs

```

# need a list of the CW and HD individuals 
fgrep -w 'CW' HD_CW_extended_list.txt | cut -f1 > CWpops
fgrep -w 'HD' HD_CW_extended_list.txt | cut -f1 > HDpops

grep -f CWpops HDCWsampleNames > CWinds.txt
grep -f HDpops HDCWsampleNames > HDinds.txt

zcat VCFs/shortName_filtered_CW_HD_SNPs.vcf.gz | vcftools --vcf - --weir-fst-pop HDinds.txt --weir-fst-pop CWinds.txt --out perSiteFST/HD_vs_CW_perSite_Fst

# this gives me a calculation of FST for all 15,299,959 SNPs 

# filtering out zero fst
fgrep -v "-" perSiteFST/HD_vs_CW_perSite_Fst.weir.fst | egrep -v -w '$0' > perSiteFST/filtered_HD_vs_CW_perSite_Fst.weir.fst

# this cuts it down to 4,569,451

```

#### Step 11. Getting per base nucleotide diversities

this was done with the no filter GL outputs
``````
realSFS saf2theta test03_allsites.saf.idx -sfs test03_allsites.sfs -outname test03_pernuc
thetaStat print test03_pernuc.thetas.idx > CW_pernuc.txt 2> /dev/null # to enable grepping specific chromosomes
``````
saf file is the output of step 6, and sfs is the output of the first part of step 7 \
after visualising this I realised that this is horrible to look at -> thus I used a sliding window approach instead: \
window length of 100 sites, and step size of 25 sites 
``````
thetaStat do_stat out.thetas.idx -win 100 -step 25  -outnames theta.thetasWindow.gz
``````
out.thetas.idx was the output of the per-site analysis, "saf2theta"

#### Step 12. Check for functional enrichment for high Fst SNPs

#### Step 13. Testing for selection at the per-SNP level using PCAngsd

(getting files: scp cpetak@pbio381.uvm.edu:/data/project_data/GroupProjects/UTR/EG_final_project/test_chr.txt .)

------    
