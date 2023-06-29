#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -J LAD_epic2
#SBATCH -p all
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%A-%a.out
#SBATCH --array=0-5
#
date
module load python3
module load ucsctools
  samples=(LMNB1_primed_r1 \
          LMNB1_primed_r2 \
          LMNB1_primed_r3 \
          LMNB1_naive_r1 \
          LMNB1_naive_r2 \
          LMNB1_naive_r3)
  #
  treatments=(s2_S2 \
    s5_S5 \
    s8_S8 \
    s11_S10 \
    s14_S13 \
    s17_S16)
  #
  controls=(S1_S1 \
    s4_S4 \
    s7_S7 \
    s10_S9 \
    s13_S12 \
    s16_S15)
  #
  bam_dir=../../LAD_analysis/naive_pipeline/LAD_bam_files
  treatment_fn="${bam_dir}/sorted_${treatments[${SLURM_ARRAY_TASK_ID}]}.bam"
  control_fn="${bam_dir}/sorted_${controls[${SLURM_ARRAY_TASK_ID}]}.bam"
  sample="${samples[${SLURM_ARRAY_TASK_ID}]}"

  mkdir results
  bin_size=10000
  gap=10
  bed_fn="results/${sample}_${bin_size}_${gap}.LAD.bed"
  epic2 --guess-bampe --genome hg38 -t "$treatment_fn" -c "$control_fn" -o "$bed_fn" \
        --bin-size "$bin_size" \
        --gaps-allowed "$gap"

  #python3 analyze.py
  tmp_fn="results/${sample}.tmp.bed"
  awk 'FNR > 1 { printf ("%s\t%s\t%s\n", $1, $2, $3)}' "$bed_fn" > "$tmp_fn"
  sorted_fn="results/${sample}.sorted.LAD.bed"
  LC_COLLATE=C sort -k1,1 -k2,2n "$tmp_fn" > "$sorted_fn" 
  bedToBigBed "$sorted_fn" results/hg38.chrom.sizes "${bed_fn}.bb"
#
module load ucsctools
module unload python3
module unload python3
module unload ucsctools
date
