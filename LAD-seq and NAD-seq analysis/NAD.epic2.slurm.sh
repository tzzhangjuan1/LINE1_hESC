#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH -J NAD_epic2
#SBATCH -p all
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%A-%a.out
#SBATCH --array=0-5
#
date
module load python3
module load ucsctools
  samples=(NAD_naive_r1 \
          NAD_naive_r2 \
          NAD_naive_r3 \
          NAD_primed_r1 \
          NAD_primed_r2 \
          NAD_primed_r3)

  #
  controls=(dedup_s1 \
    dedup_s3 \
    dedup_s5 \
    dedup_s7 \
    dedup_s9 \
    dedup_s11)
  #
  treatments=(dedup_s2 \
    dedup_s4 \
    dedup_s6 \
    dedup_s8 \
    dedup_s10 \
    dedup_s12)
  #
  bam_dir=../../NAD_bam_files
  treatment_fn="${bam_dir}/${treatments[${SLURM_ARRAY_TASK_ID}]}.bam"
  control_fn="${bam_dir}/${controls[${SLURM_ARRAY_TASK_ID}]}.bam"
  sample="${samples[${SLURM_ARRAY_TASK_ID}]}"

  mkdir results
  bin_size=10000
  gap=5
  bed_fn="results/${sample}_${bin_size}_${gap}.NAD.bed"
  #epic2 --guess-bampe --genome hg38 -t "$treatment_fn" -c "$control_fn" -o "$bed_fn" \
  #      --bin-size "$bin_size" \
  #      --gaps-allowed "$gap"
  bw_fn="results/${sample}_${bin_size}_${gap}.bw"
  epic2-bw --version --genome hg38 -t "$treatment_fn" -c "$control_fn" \
        --bin-size "$bin_size" \
        --log2fc-bigwig "$bw_fn"
  #
  #tmp_fn="results/${sample}.tmp.bed"
  #awk 'FNR > 1 { printf ("%s\t%s\t%s\n", $1, $2, $3)}' "$bed_fn" > "$tmp_fn"
  #sorted_fn="results/${sample}.sorted.NAD.bed"
  #LC_COLLATE=C sort -k1,1 -k2,2n "$tmp_fn" > "$sorted_fn" 
  #bedToBigBed "$sorted_fn" results/hg38.chrom.sizes "${bed_fn}.bb"
#
module unload python3
module unload ucsctools
date
