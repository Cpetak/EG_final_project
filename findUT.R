####### Find UTRs based on depth of coverage of 3' RNA-seq reads
##### Input file is output of samtools depth command using bed file of reference seqences with known coding domains removed 
### Output is a bed file of locations of UTRs
## BWW
# 4/24/20


library(plyr)
library(tidyverse)
library(magrittr)

depthDat <- read_csv("groupProj/final_RNAbamDepth.csv", col_names = c("contig", "base", "depth"))

# this file has already had sites with depth of 0-3 filtered out, but might need to filter more

depthDat %>%
  filter(depth < 5) %>%
  nrow() # 2127646  13% of the sites

# just leave it at 5 for now 

depthDat %<>%
  filter(depth >= 5)

test <- filter(depthDat, contig=="MA_100009")

# need helper function

seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 


# another helper function to pull out sequences from seql output
# where x is one index value of the seql output
procSeq <- function(x,seqs,name){
  
  return(data.frame(contig=name,
                    start = seqs$values[x],  end=seqs$values[x]+(seqs$lengths[x]-1)))
  
}



allUTR <- function(x) {
  
  name <- as.character(x[1,1])
  # find consecutive stretechs of coverage within the contig 
  seqs <- seqle(x$base)
  
  # get the index for all the stretches of coverage over 150bp
  utrKey <- which(seqs$lengths  >150)
  
  # now build the output df made up of the contig and pos for each strech over 150
  utrKey %>%
    map_df(procSeq, seqs=seqs, name=name) # see above
  
  # this should return a DF with one row for each stretch of coverage over 150bp
}


# test with data subset (single contig)
allUTR(test)

# apply across all contigs

utrFinal <- ddply(depthDat, "contig", allUTR, .drop = T)

nrow(utrFinal) #22220 which is not too far under the 28,354 HC gene predictions

UTRsizesFinal <- utrFinal %>% 
  mutate(size= end - start) %>%
  select(contig, size)

quantile(UTRsizesFinal$size)
#0%  25%  50%  75% 100% 
#150  181  235  350 7493 

# make bed file 

utrFinal$label <- rep("UTR", nrow(utrFinal))

# this needs to be zero-based for annotation 

utrFinal %<>%
  mutate(start=start-1)

write_tsv(utrFinal, "finalUTRs.bed", col_names = F)
