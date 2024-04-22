#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH -N 1
#SBATCH -n 6


# 1.	Sort bed file (otherwise cannot intersect)
bedtools sort -i 2312simple_naive_LAD.bed > 2312simple_naive_LAD_sorted.bed
bedtools sort -i 2312simple_naive_NAD.bed > 2312simple_naive_NAD_sorted.bed
bedtools sort -i 2312simple_primed_LAD.bed > 2312simple_primed_LAD_sorted.bed
bedtools sort -i 2312simple_primed_NAD.bed > 2312simple_primed_NAD_sorted.bed


# 2.	Intersect bed file (to co-LADsNADs bed)
 bedtools intersect -a 2312simple_naive_LAD_sorted.bed -b 2312simple_naive_NAD_sorted.bed > 2312naive_LAD-NAD_Overlaps.bed
 bedtools intersect -a 2312simple_primed_LAD_sorted.bed -b 2312simple_primed_NAD_sorted.bed > 2312primed_LAD-NAD_Overlaps.bed

# 3.	Multiinter (Bed covers both -A.bed and -B.bed)
 bedtools multiinter -header -i 2312simple_naive_NAD_sorted.bed 2312simple_naive_LAD_sorted.bed > 2312simple_naiveLAD-NAD_multiinter.bed
 bedtools multiinter -header -i 2312simple_primed_LAD_sorted.bed 2312simple_primed_NAD_sorted.bed > 2312simple_primedLAD-NAD_multiinter.bed

# 4.	Complement 
 bedtools sort -i hg38_simple.genome> sorted_hg38_simple.genome

 bedtools complement -i 2312simple_naiveLAD-NAD_multiinter2.bed -g sorted_hg38_simple.genome > 2312naive_nonLAD-NAD.bed
 bedtools complement -i 2312simple_primedLAD-NAD_multiinter2.bed -g sorted_hg38_simple.genome > 2312primed_nonLAD-NAD.bed


#5. subtract (to get LADs-only, NADs-only bed)
bedtools subtract -a 2312simple_naive_LAD_sorted.bed -b 2312naive_LAD-NAD_Overlaps.bed > 2312naive_LAD-only.bed
bedtools subtract -a 2312simple_naive_NAD_sorted.bed -b 2312naive_LAD-NAD_Overlaps.bed > 2312naive_NAD-only.bed
bedtools subtract -a 2312simple_primed_LAD_sorted.bed -b 2312primed_LAD-NAD_Overlaps.bed > 2312primed_LAD-only.bed
bedtools subtract -a 2312simple_primed_NAD_sorted.bed -b 2312primed_LAD-NAD_Overlaps.bed > 2312primed_NAD-only.bed

 
# 6. intersect for LADs,NADs bed overlapped genes
 bedtools intersect -a 2312simple_naive_LAD_sorted.bed -b hg38_genes_refseq.bed -wa -wb > 2312naive_LAD_hg38ref.bed             
 bedtools intersect -a 2312simple_naive_NAD_sorted.bed -b hg38_genes_refseq.bed -wa -wb > 2312naive_NAD_hg38ref.bed
 bedtools intersect -a 2312simple_primed_LAD_sorted.bed -b hg38_genes_refseq.bed -wa -wb > 2312primed_LAD_hg38ref.bed
 bedtools intersect -a 2312simple_primed_NAD_sorted.bed -b hg38_genes_refseq.bed -wa -wb > 2312primed_NAD_hg38ref.bed

 ## edit in txt to get genes list, and remove duplicated genes
 ## LAD-only, NAD-only, co-LADNAD and non-LADNAD genes obtained in R. 



done

exit 0




