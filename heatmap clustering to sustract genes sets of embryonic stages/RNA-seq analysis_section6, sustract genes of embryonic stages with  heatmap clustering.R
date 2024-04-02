# This code is to input the normalized counts table form the scRNA-seq by Yan et al., 2013, to subtract clustered genes set at specific embryonic stages. 
# Genes of each embryos stages are listed Table S2.  

# The data resource from Yan et al. are available from their Supplemental Table.
# The input "scRNAseq_Yan, human embryo stages.txt" tables above is available in the repository. 


library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)

############################## heat map plot ###################################

topT <- read.table("scRNAseq_Yan, human embryo stages.txt", header=T, sep="\t",comment.char="#")
subtopT <- subset(topT, rowSums(topT[,2:12])>10)   # remove sum-each-row <10
rnames <- subtopT[,1]
mat_subtopT <- data.matrix(subtopT[,2:ncol(topT)])
rownames(mat_subtopT) <- rnames
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),               # for blue
               seq(0.01,0.8,length=100),           # for white
               seq(0.81,1,length=100))             # for red
## Row clustering (adjust here distance/linkage methods to what you need!)
hr1 <- hclust(as.dist(1-cor(t(mat_subtopT), method="pearson")), method="complete")
#hr2 <- hclust(as.dist(1-cor(t(mat_subtopT), method="spearman")), method="complete")
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(as.dist(1-cor(mat_subtopT, method="spearman")), method="complete")
png("../toptable_heatmap1.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*900,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(mat_subtopT,
          #cellnote = mat_data,         # same data set for cell labels
          main = "human embryo stages,pearson",   # heat map title
          #notecol="black",             # change font color of cell labels to black
          notecol="NA",
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          scale="row",
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="col",     # only draw a col dendrogram
          Rowv = as.dendrogram(hr1), # apply default clustering method
          Colv=as.dendrogram(hc))           # turn off column clustering

# the heatmap output displays clustered genes expressed at specific embryonic stages)


############################## extract dendrogram cluster from pheatmap ###################################

out = heatmap.2(mat_subtopT,
                #cellnote = mat_data,         # same data set for cell labels
                main = "human embryo stages,pearson",   # heat map title
                #notecol="black",             # change font color of cell labels to black
                notecol="NA",
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(12,9),     # widens margins around plot
                col=my_palette,       # use on color palette defined earlier
                scale="row",
                #breaks=col_breaks,    # enable color transition at specified limits
                dendrogram="col",     # only draw a col dendrogram
                Rowv = as.dendrogram(hr1), # apply default clustering method
                Colv=as.dendrogram(hc))

################################ cut-tree on column ###################################
clus_col <- cutree(as.hclust(out$colDendrogram), 3)
names(clus_col[clus_col==1])
# [1] "X8cell"  "Morulae"
names(clus_col[clus_col==2])
# [1] "Oocyte" "Zygote" "X2cell" "X4cell"
names(clus_col[clus_col==3])
# [1] "Trophectoderm"  "Primitive"      "epiblast"       "hESC.p0"        "hESC.passage10"

### In this way it subtract sample names, by "cut the column tree"

### Similar as above, if "cut the row clustering tree", genes of the same cluster are subtraced
################################ cut-tree on row ###################################
# plot the raw tree (genes cluster)
plot(hr1)   
# this is the raw (genes) cluster in Extended Data Fig.5m. 

##### use abline() to find the "cut tree site", so subtract clusters of interest ##### 
plot(hr1)
abline(h = 1.92, col = "grey30", lty = 2, lwd = 2)
abline(h = 1.95, col = "red2", lty = 2, lwd = 2)
abline(h = 1.86, col = "blue2", lty = 2, lwd = 2)

## cut-tree on row based on the abline() above, and subtract genes in each cluster #####
clus_row <- sort(cutree(as.hclust(out$rowDendrogram),h = 1.86))
names(clus_row[clus_row==1])     # the genes and genes number (2457) in cluster 1 
names(clus_row[clus_row==2])     # the genes and genes number (3930) in cluster 2 
names(clus_row[clus_row==3])     # the genes and genes number (433) in cluster 3
names(clus_row[clus_row==4])     # the genes and genes number (1687) in cluster 4
names(clus_row[clus_row==5])     # the genes and genes number (546) in cluster 5
names(clus_row[clus_row==6])     # the genes and genes number (2721) in cluster 6
names(clus_row[clus_row==7])     # the genes and genes number (211) in cluster 7
names(clus_row[clus_row==8])     # the genes and genes number (301) in cluster 8
names(clus_row[clus_row==9])     # the genes and genes number (291) in cluster 9
# Based on the cluster size/position, try to figure out which cluster mach to which embryonic stage.

# export these clustered genes to a csv file, in the same order as showing in the dendrogram,the order is bottom-to-top 
write.table(data.frame(gene = names(clus_row[clus_row==1])),'cluster1.csv',row.names = FALSE,quote = FALSE,sep = ',') #export cluster 1 that contains 2457 genes
# Export genes in the other clusters, and see Extended Data Fig.5m for the corresponding embryonic stage. 

# Genes sets extracted by this analysis are listed in Supplementary Table 2. 


