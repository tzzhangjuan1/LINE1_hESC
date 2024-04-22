#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=40G
#SBATCH -J NAD_bam
#SBATCH -p all
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%A-%a.out
#SBATCH --array=1-12
#
date
module load ucsctools
module load deeptools

#peak_fn="../../NAD_finder_R/s${SLURM_ARRAY_TASK_ID}.called.broad.peak.bed"
#sorted_peak_fn="results/s${SLURM_ARRAY_TASK_ID}.sort.bed"
#if [[ -f "$peak_fn" ]]; then
#  LC_COLLATE=C sort -k1,1 -k2,2n "${peak_fn}" > ${sorted_peak_fn}
#  bedToBigBed "${sorted_peak_fn}" results/hg38.chrom.sizes "results/s${SLURM_ARRAY_TASK_ID}.bb"
#fi
#bin_size=10000
#bam_fn="../../NAD_bam_files/dedup_s${SLURM_ARRAY_TASK_ID}.bam"
#bamCoverage --bam ${bam_fn} -o "results/s${SLURM_ARRAY_TASK_ID}_${bin_size}.bw" \
#    --binSize ${bin_size} \
#    --normalizeUsing RPGC \
#    --effectiveGenomeSize 2913022398 \
#    --extendReads
####    --ignoreForNormalization chrX \
#
#r=$((${SLURM_ARRAY_TASK_ID}%2))
#if [ $r -gt 0 ]; then
#  bin_size=10000
#  n=$((${SLURM_ARRAY_TASK_ID}+1))
#  bw_fn_1="results/s${n}_${bin_size}.bw"
#  bw_fn_2="results/s${SLURM_ARRAY_TASK_ID}_10000.bw"
#  bigwigCompare --bigwig1 ${bw_fn_1} --bigwig2 ${bw_fn_2} -o "results/log2_s${SLURM_ARRAY_TASK_ID}_${bin_size}.bw"
#fi
#

#####
# 1bp


module unload ucsctools
module unload deeptools
date
