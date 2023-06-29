#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=80G
#SBATCH -J NAD_bam
#SBATCH -p all
#SBATCH -c 16
#SBATCH -N 1
#SBATCH -o %x-%A-%a.out
#SBATCH --array=4-4
#
date
module load bwa
module load samtools/1.10
module load picard
ref_file=../RepeatHiC/data/2020_02_19/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#
fn="s${SLURM_ARRAY_TASK_ID}"
ori_bam_fn="NAD_bam_files/ori_${fn}.bam"
# -S: ignore incompatibility with previous version
# -F 4: Unmapped read
# -q 20: MAPQ >= 20
# -M: Mark shorter split hits as secondary (for Picard compatibility).
bwa mem -t 16 -M ${ref_file} "NAD_data_split/${fn}.R1.fastq.gz" "NAD_data_split/${fn}.R2.fastq.gz" \
| samtools view -b -S -o "$ori_bam_fn" -
#
filtered_bam_fn="NAD_bam_files/filtered_${fn}.bam"
samtools view -b -F 4 -q 20 "$ori_bam_fn" > "$filtered_bam_fn" 
#
sorted_fn="NAD_bam_files/sorted_${fn}.bam"
samtools sort "$filtered_bam_fn" > "$sorted_fn" 
#dedup
dup_fn="NAD_bam_files/dup_${fn}.bam"
metric_fn="NAD_bam_files/met_${fn}.txt"
java -jar ${picard_dir}/picard.jar MarkDuplicates I="$sorted_fn" O="$dup_fn" M="$metric_fn"
samtools index "$dup_fn"
#Expunge marked duplicate reads, and then index new BAM
dedup_fn="NAD_bam_files/dedup_${fn}.bam"    
samtools view -b -F 0x400 "$dup_fn" > "$dedup_fn"
samtools index "$dedup_fn"
#
module unload bwa
module unload samtools
module unload picard
date
