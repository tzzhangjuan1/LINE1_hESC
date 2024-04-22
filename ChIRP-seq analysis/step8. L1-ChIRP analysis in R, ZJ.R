library(limma)
library(gdata)
library(gplots)
library(tidyverse)
library(ggrepel) 
library(viridis)
library(ggpubr)

library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

########## calculate L1ChIRP length, by chromosomes ###########
L1ChIRP <- as.data.frame(read.table("broad_l1-chirp_merge_peaks_blacklisted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

L1ChIRP_Length <- L1ChIRP %>% 
  group_by(V1) %>%
  mutate(length=V3-V2) %>%
  summarise(Sum=sum(length))
write.csv(L1ChIRP_Length,"L1ChIRP_Length.csv", row.names = FALSE)

########## calculate L1ChIRP length ###########

L1HSPA_Overlaps <- as.data.frame(read.table("broad_L1-chirp_L1HSPA-overlap.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
L1HSPA_w0to500 <- as.data.frame(read.table("broad_L1-chirp_L1HSPA-w0to500.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))  
L1HSPA_w500to2000 <- as.data.frame(read.table("broad_L1-chirp_L1HSPA-w500to2000.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
L1HSPA_out_w2000 <- as.data.frame(read.table("broad_L1-chirp_L1HSPA-out_w2000.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))  


sum(L1HSPA_Overlaps$V3-L1HSPA_Overlaps$V2)         # L1HSPA_Overlaps_Length total, 110303222
sum(L1HSPA_w0to500$V3-L1HSPA_w0to500$V2)         # L1HSPA_w0to500_Length total, 247040978
sum(L1HSPA_w500to2000$V3-L1HSPA_w500to2000$V2)   # L1HSPA_w500to2000_Length total,  6010417
sum(L1HSPA_out_w2000$V3-L1HSPA_out_w2000$V2)     # L1HSPA_out_w2000, 22652945

########## calculate L1ChIRP-genes numbers in LADs/NADs ###########
L1ChIRP_gene <- read.table("L1ChIRP-genes.txt",header = FALSE, sep="\t")    #5851

naive_LAD_only <- read.table("2312naive_LAD_only_hg38ref_genes_rmDup_R.txt", header=F, sep="\t")      #1686
LAD_only_L1ChIRP <- naive_LAD_only[naive_LAD_only$V1 %in% L1ChIRP_gene$V1,]         #317

naive_NAD_only <- read.table("2312naive_NAD_only_hg38ref_genes_rmDup_R.txt", header=F, sep="\t")      #14499
NAD_only_L1ChIRP <- naive_NAD_only[naive_NAD_only$V1 %in% L1ChIRP_gene$V1,]         #2323

naive_coLADNAD <- read.table("2312naive_coLADNAD_hg38ref_genes_rmDup_R.txt", header=F, sep="\t")      #5637
coLADNAD_L1ChIRP <- naive_coLADNAD[naive_coLADNAD$V1 %in% L1ChIRP_gene$V1,]         #2914

naive_nonLADNAD <- read.table("2312naive_nonLADNAD_hg38ref_genes_rmDup_R.txt", header=F, sep="\t")      #6489
nonLADNAD_L1ChIRP <- naive_nonLADNAD[naive_nonLADNAD$V1 %in% L1ChIRP_gene$V1,]         #297


