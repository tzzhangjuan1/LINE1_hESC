#!/bin/bash
#SBATCH --job-name=bedwindow

# 1.	Sort bed file (otherwise cannot intersect)
# sort -k1,1 -k2,2n broad_l1-chirp_merge_peaks_blacklisted.bed > broad_l1-chirp_merge_peaks_blacklisted_sort.bed
# sort -k1,1 -k2,2n hg38_2020_rmsk_onlyHumanL1.bed > L1HSL1PAs_sort.bed

# 2.	Intersect bed file (B region overlap with A)
# bedtools intersect -a L1HSL1PAs_sort.bed -b broad_l1-chirp_merge_peaks_blacklisted_sort.bed > broad_L1-chirp_L1HSPA-overlap.bed

# 2.	Intersect bed file (B region within -w of A)
bedtools window -w 500 -a L1HSL1PAs_sort.bed -b broad_l1-chirp_merge_peaks_blacklisted_sort.bed > broad_L1-chirp_L1HSPA-w500.bed     #open as txt, rm L1HAPA loci
bedtools window -w 1000 -a L1HSL1PAs_sort.bed -b broad_l1-chirp_merge_peaks_blacklisted_sort.bed > broad_L1-chirp_L1HSPA-w1000.bed   #open as txt, rm L1HAPA loci
bedtools window -w 2000 -a L1HSL1PAs_sort.bed -b broad_l1-chirp_merge_peaks_blacklisted_sort.bed > broad_L1-chirp_L1HSPA-w2000.bed   #open as txt, rm L1HAPA loci

#5. subtract (to get 0-500, 500-1000, 1000-2000 bed)
# bedtools sort -i broad_L1-chirp_L1HSPA-w500_rmL1.bed > broad_L1-chirp_L1HSPA-w500_rmL1_sort.bed
# bedtools sort -i broad_L1-chirp_L1HSPA-w1000_rmL1.bed > broad_L1-chirp_L1HSPA-w1000_rmL1_sort.bed
# bedtools sort -i broad_L1-chirp_L1HSPA-w2000_rmL1.bed> broad_L1-chirp_L1HSPA-w2000_rmL1_sort.bed

# bedtools subtract -a broad_L1-chirp_L1HSPA-w500_rmL1.bed -b broad_L1-chirp_L1HSPA-overlap.bed > broad_L1-chirp_L1HSPA-w0to500.bed

# bedtools subtract -a broad_L1-chirp_L1HSPA-w1000_rmL1.bed -b broad_L1-chirp_L1HSPA-w500_rmL1.bed > broad_L1-chirp_L1HSPA-w500to1000.bed
# bedtools subtract -a broad_L1-chirp_L1HSPA-w2000_rmL1.bed -b broad_L1-chirp_L1HSPA-w1000_rmL1.bed > broad_L1-chirp_L1HSPA-w1000to2000.bed

# bedtools subtract -a broad_L1-chirp_L1HSPA-w2000_rmL1.bed -b broad_L1-chirp_L1HSPA-w500_rmL1.bed > broad_L1-chirp_L1HSPA-w500to2000.bed
# bedtools subtract -a broad_l1-chirp_merge_peaks_blacklisted_sort.bed -b broad_L1-chirp_L1HSPA-w2000_rmL1.bed > broad_L1-chirp_L1HSPA-out_w2000.bed

# 4. Complement (to get >2000 bed)
# bedtools sort -i hg38_simple.genome> sorted_hg38_simple.genome
# bedtools complement -i broad_L1-chirp_L1HSPA-w2000_rmL1_sort.bed -g broad_l1-chirp_merge_peaks_blacklisted_sort.bed > broad_L1-chirp_L1HSPA-out_w2000.bed



done
exit 0
