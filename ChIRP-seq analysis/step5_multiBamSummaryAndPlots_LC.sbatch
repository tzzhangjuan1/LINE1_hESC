#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=10
#SBATCH --job-name=multiBamSummary

###these tools are just part of bedtools still
export PATH=/home/wranalab/lcaldwel/programs/bedtools:$PATH


#Example
#multiBamSummary bins \
 #--bamfiles testFiles/*bam \ # using all BAM files in the folder
 #--minMappingQuality 30 \
 #--region 19 \ # limiting the binning of the genome to chromosome 19
 #--labels H3K27me3 H3K4me1 H3K4me3 HeK9me3 input \
 #-out readCounts.npz --outRawCounts readCounts.tab

#should be able to just call these, don't need to reference bedtools specifically?
#this set of commands doesn't work with number of processors command, it just gives a bin error ????

echo "starting multibam summary of malat1-ChIRP..."
multiBamSummary bins \
--bamfiles filtered_s2.bam filtered_s4.bam \
--labels malat1-ChIRP-rep1 malat1-ChIRP-rep2 \
-o malat1-ChIRP_rep1and2.npz --outRawCounts malat1-ChIRP-readCounts.tab -p 6 #-p 6#--numberOfProcessors 8

echo "starting multibam summary of L1-ChIRP..."
multiBamSummary bins \
--bamfiles filtered_s3.bam filtered_s5.bam \
--labels L1-ChIRP-rep1 L1-ChIRP-rep2  \
-o L1-ChIRP_rep1and2.npz --outRawCounts L1-ChIRP-readCounts.tab -p 6  #-p 6#--numberOfProcessors 8

echo "starting multibam summary of K27me3-ChIRP-all..."
multiBamSummary bins \
--bamfiles filtered_S1.bam filtered_s2.bam filtered_s3.bam filtered_s4.bam filtered_s5.bam \
--labels L1-ChIRP-input malat1-ChIRP-rep1 L1-ChIRP-rep1 malat1-ChIRP-rep2 L1-ChIRP-rep2 \
-o K27me3-ChIRP-all-S1to5.npz --outRawCounts K27me3-ChIRP-all-S1to5-readCounts.tab -p 6 #
#uses npz file for PCA


echo "plotting PCA..."
#this is one that I added
plotPCA -in K27me3-ChIRP-all-S1to5.npz \
-o K27me3-ChIRP-all_S1-5_PCA_readCounts.png \
-T "PCA of read counts, K27me3-ChIRP-all S1-S5"

echo "end of PCA..."

###plotFingerprint
#This quality control will most likely be of interest for you if you are dealing with ChIRP-seq
#samples as a pressing question in ChIRP-seq experiments is “Did my ChIRP work?”, i.e. did the antibody-treatment
#enrich sufficiently so that the ChIRP signal can be separated from the background signal?
echo "plotting fingerprints..."
#(After all, around 90% of all DNA fragments in a ChIRP experiment will represent the genomic background).

plotFingerprint -b filtered_S1.bam filtered_s2.bam filtered_s3.bam filtered_s4.bam filtered_s5.bam \
--labels L1-ChIRP-input malat1-ChIRP-rep1 L1-ChIRP-rep1 malat1-ChIRP-rep2 L1-ChIRP-rep2 \
-T "Fingerprints of K27me3-ChIRP all samples S1-S5"  \
--plotFile plotFingerprint_K27me3-ChIRP-all-S1to5.png \
--outRawCounts plotFingerprint_K27me3-ChIRP-all-S1to5.tab

echo "plotting correlation scatterplot..."
#scatterplot
plotCorrelation \
-in K27me3-ChIRP-all-S1to5.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of K27me3-ChIRP-all S1-S5" \
--whatToPlot scatterplot \
-o scatterplot_K27me3-ChIRP-all_PearsonCorr_bigwigScores-S1to5.png   \
--outFileCorMatrix K27me3-ChIRP-all_PearsonCorr_bigwigScores-S1to5.tab

echo "plotting correlation heatmap..."
#heatmap
plotCorrelation \
    -in K27me3-ChIRP-all-S1to5.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of K27me3-ChIRP-all S1-S5" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o K27me3-ChIRP-all_heatmap_SpearmanCorr_readCounts-S1-5.png   \
    --outFileCorMatrix K27me3-ChIRP-all_heatmap_SpearmanCorr_readCounts-S1-5.tab



echo "Finished ChIRPseq step 5: multi-bam summary and plots with exit code $? at: `date`"
