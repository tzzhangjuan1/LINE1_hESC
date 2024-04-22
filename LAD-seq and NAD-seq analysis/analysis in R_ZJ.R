
library(limma)
library(gdata)
library(gplots)
library(tidyverse)

library(edgeR)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ggpubr)

library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

setwd("~/")   #run code in the same folder.

#####################################################
##### code Run in R, to get a clear genes list  ##### 
#####################################################

hg38_Gene <- read.table("UCSC RefSeq_hg38_rmDup.txt", header=F, sep="\t")      #28307 
naive_LAD <- read.table("2312naive_LAD_hg38ref_genes_rmDup.txt", header=F, sep="\t")      #7321
naive_NAD <- read.table("2312naive_NAD_hg38ref_genes_rmDup.txt", header=F, sep="\t")      #20134
naive_coLADNAD <- subset(naive_LAD, naive_LAD$V1 %in% naive_NAD$V1)     #5636
naive_LAD_only <- subset(naive_LAD, !(naive_LAD$V1 %in% naive_coLADNAD$V1))    #1685
naive_NAD_only <- subset(naive_NAD, !(naive_NAD$V1 %in% naive_coLADNAD$V1))    #14498
naive_nonLADNAD <- subset(hg38_Gene, !(hg38_Gene$V1 %in% naive_NAD$V1) & !(hg38_Gene$V1 %in% naive_LAD$V1))  #6488

write.csv(naive_LAD_only,"2312naive_LAD_only_hg38ref_genes_rmDup_R.txt", row.names = FALSE)
write.csv(naive_NAD_only,"2312naive_NAD_only_hg38ref_genes_rmDup_R.txt", row.names = FALSE)
write.csv(naive_coLADNAD,"2312naive_coLADNAD_hg38ref_genes_rmDup_R.txt", row.names = FALSE)
write.csv(naive_nonLADNAD,"2312naive_nonLADNAD_hg38ref_genes_rmDup_R.txt", row.names = FALSE)

#######################################################
######## calculate LADsNADs_coverage _length ##########
#######################################################

# read in "2312naive_LAD-NAD_Overlaps.bed" for example:
Naive_LAD_NAD_Overlaps <- as.data.frame(read.table("2312naive_LAD-NAD_Overlaps.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))  

sum(Naive_LAD_NAD_Overlaps$V3-Naive_LAD_NAD_Overlaps$V2)         # Naive_LAD_NAD_Overlaps total length 


## calculate each chromosome length

Naive_LAD_NAD_Overlaps_Length <- Naive_LAD_NAD_Overlaps %>% 
  group_by(V1) %>%
  mutate(length=V3-V2) %>%
  summarise(Sum=sum(length))
write.csv(Naive_LAD_NAD_Overlaps_Length,"Naive_LAD_NAD_Overlaps.csv", row.names = FALSE)

