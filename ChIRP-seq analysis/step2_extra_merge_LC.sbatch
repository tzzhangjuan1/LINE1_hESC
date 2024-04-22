#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name chip1.5

##bowtie2 cant run in parallel so don't request more nodes
export PATH=/home/wranalab/lcaldwel/programs/samtools-1.18:$PATH
export PATH=/home/wranalab/lcaldwel/.local/bin/deeptools:$PATH


###merge with samtools
echo "merging with samtools..."
#this section failed the first time because the files came out of part 1 as _tr_s2.bam etc
#samtools merge malat1-chirp_merge.bam filtered_s2.bam filtered_s4.bam
#samtools merge l1-chirp_merge.bam filtered_s3.bam filtered_s5.bam
samtools index malat1-chirp_merge.bam
#samtools index l1-chirp_merge.bam

#get coverage
  #this command is from DEEP TOOLS
bamCoverage -b malat1-chirp_merge.bam -o malat1-chirp_merge.bw -v --normalizeUsing CPM
bamCoverage -b l1-chirp_merge.bam -o l1-chirp_merge.bw -v --normalizeUsing CPM

echo "Finished chipseq step 3: Merging and deeptools run with exit code $? at: `date`"
