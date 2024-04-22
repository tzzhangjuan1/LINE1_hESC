#!/bin/bash
#SBATCH --job-name=bamCoverage
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB

export PATH=/home/wranalab/lcaldwel/programs/bedtools:$PATH
#uses deeptools which I installed with pip
blacklist=/home/wranalab/lcaldwel/ChipSeq2024/blacklist/hg38.Kundaje.GRCh38_unified_blacklist.bed

mkdir graphicouts_step5

echo "Starting bam coverage..."

###get bam coverage (formerly bamCoverage.sh)
bamCoverage -b filtered_S1.bam -o ctr-input_RPKMnorm.bw -v --normalizeUsing RPKM
bamCoverage -b filtered_s2.bam -o MALAT1-ChIRP-rep1-RPKMnorm.bw -v --normalizeUsing RPKM
bamCoverage -b filtered_s3.bam -o L1-ChIRP-rep1-RPKMnorm.bw -v --normalizeUsing RPKM
bamCoverage -b filtered_s4.bam -o MALAT1-ChIRP-rep2-RPKMnorm.bw -v --normalizeUsing RPKM
bamCoverage -b filtered_s5.bam -o L1-ChIRP-rep2-RPKMnorm.bw -v --normalizeUsing RPKM

##merged
bamCoverage -b malat1-chirp_merge.bam -o MALAT1-ChIRP_merge_RPKMnorm.bw -v --normalizeUsing RPKM
bamCoverage -b l1-chirp_merge.bam -o L1-ChIRP_merge_RPKMnorm.bw -v --normalizeUsing RPKM


echo "Starting bigwig averages..."
###get averages (formerly bigwigAverage.sh)

#this is all non control samples per set
bigwigAverage -b  MALAT1-ChIRP-rep1-RPKMnorm.bw MALAT1-ChIRP-rep2-RPKMnorm.bw -o Average_MALAT1-ChIRP-rep12.bw
bigwigAverage -b L1-ChIRP-rep1-RPKMnorm.bw L1-ChIRP-rep2-RPKMnorm.bw -o Average_L1-ChIRP-rep12.bw


echo "Finished chipseq step 6: normalize with RPKM, bigwig average, + bigwig compare + bam coverage bigwigs with exit code $? at: `date`"
