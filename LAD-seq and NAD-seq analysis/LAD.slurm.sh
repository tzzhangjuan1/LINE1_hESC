#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=80G
#SBATCH -J LAD_pipeline
#SBATCH -p all
#SBATCH -c 8
#SBATCH -N 1
#SBATCH -o %x-%A-%a.out
#SBATCH --array=0-9
#
date
module load bwa
module load samtools/1.10
module load ucsctools
module load deeptools

ref_file=/mnt/work1/users/hoffmangroup/lhuynh/RepeatHiC/data/2020_02_19/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#
bam_dir=LAD_bam_files
mkdir "${bam_dir}"
fastq_dir=../../LAD_raw_data/HX5EH4X/ZHA17556.20210618/210617_D00430_0462_AHLNV7BCX3/
fastqs=(S1_S1_L001_R1_001.fastq.gz \
 S1_S1_L002_R1_001.fastq.gz \
 s2_S2_L001_R1_001.fastq.gz \
 s2_S2_L002_R1_001.fastq.gz \
 s3_S3_L001_R1_001.fastq.gz \
 s3_S3_L002_R1_001.fastq.gz \
 s4_S4_L001_R1_001.fastq.gz \
 s4_S4_L002_R1_001.fastq.gz \
 s5_S5_L001_R1_001.fastq.gz \
 s5_S5_L002_R1_001.fastq.gz \
 s6_S6_L001_R1_001.fastq.gz \
 s6_S6_L002_R1_001.fastq.gz \
 s7_S7_L001_R1_001.fastq.gz \
 s7_S7_L002_R1_001.fastq.gz \
 s8_S8_L001_R1_001.fastq.gz \
 s8_S8_L002_R1_001.fastq.gz \
 s10_S9_L001_R1_001.fastq.gz \
 s10_S9_L002_R1_001.fastq.gz \
 s11_S10_L001_R1_001.fastq.gz \
 s11_S10_L002_R1_001.fastq.gz \
 s12_S11_L001_R1_001.fastq.gz \
 s12_S11_L002_R1_001.fastq.gz \
 s13_S12_L001_R1_001.fastq.gz \
 s13_S12_L002_R1_001.fastq.gz \
 s14_S13_L001_R1_001.fastq.gz \
 s14_S13_L002_R1_001.fastq.gz \
 s15_S14_L001_R1_001.fastq.gz \
 s15_S14_L002_R1_001.fastq.gz \
 s16_S15_L001_R1_001.fastq.gz \
 s16_S15_L002_R1_001.fastq.gz \
 s17_S16_L001_R1_001.fastq.gz \
 s17_S16_L002_R1_001.fastq.gz)
samples=(S1_S1 \
 s2_S2 \
 s3_S3 \
 s4_S4 \
 s5_S5 \
 s6_S6 \
 s7_S7 \
 s8_S8 \
 s10_S9 \
 s11_S10 \
 s12_S11 \
 s13_S12 \
 s14_S13 \
 s15_S14 \
 s16_S15 \
 s17_S16)

fn="${fastqs[${SLURM_ARRAY_TASK_ID}]}"
echo "$fn"
ori_bam_fn="${bam_dir}/ori_${fn}.bam"
# -S: ignore incompatibility with previous version
# -F 4: Unmapped read
# -q 10: MAPQ >= 10
# -M: Mark shorter split hits as secondary (for Picard compatibility).
#
########## Alignment
#bwa mem -t 8 -M ${ref_file} "${fastq_dir}/${fn}"\
#| samtools view -b -S -o "$ori_bam_fn" -
#
filtered_bam_fn="${bam_dir}/filtered_${fn}.bam"
#samtools view -b -F 4 -q 10 "$ori_bam_fn" > "$filtered_bam_fn" 
#
sorted_fn="${bam_dir}/sorted_${fn}.bam"
#samtools sort "$filtered_bam_fn" > "$sorted_fn"
#samtools index "$sorted_fn"
#
########## Merge replicates 
sample="${samples[${SLURM_ARRAY_TASK_ID}]}"
merged_fn="${bam_dir}/sorted_${sample}.bam"
#samtools merge "${merged_fn}" \
#  "${bam_dir}/sorted_${sample}_L001_R1_001.fastq.gz.bam" \
#  "${bam_dir}/sorted_${sample}_L002_R1_001.fastq.gz.bam" 
#samtools index "${merged_fn}"
########## Export to bigwig files
bin_size=10000
#bamCoverage --bam "${merged_fn}" -o "${bam_dir}/${sample}_${bin_size}.bw" \
#    --binSize "${bin_size}" \
#    --normalizeUsing RPGC \
#    --effectiveGenomeSize 2913022398
#bamCoverage --bam "${merged_fn}" -o "${bam_dir}/${sample}_CPM_${bin_size}.bw" \
#    --binSize "${bin_size}" \
#    --normalizeUsing CPM \
#    --effectiveGenomeSize 2913022398

########## Export to log fold change
treatments=(s2_S2 \
 s3_S3 \
 s5_S5 \
 s6_S6 \
 s8_S8 \
 s11_S10 \
 s12_S11 \
 s14_S13 \
 s15_S14 \
 s17_S16)

controls=(S1_S1 \
 S1_S1 \
 s4_S4 \
 s4_S4 \
 s7_S7 \
 s10_S9 \
 s10_S9 \
 s13_S12 \
 s13_S12 \
 s16_S15)

treatment="${treatments[${SLURM_ARRAY_TASK_ID}]}"
control="${controls[${SLURM_ARRAY_TASK_ID}]}"

#bigwigCompare --bigwig1 ${bam_dir}/${treatment}_${bin_size}.bw --bigwig2 ${bam_dir}/${control}_${bin_size}.bw \
#  --operation ratio -o "${bam_dir}/ratio_${treatment}_${bin_size}.bw"

bigwigCompare --bigwig1 ${bam_dir}/${treatment}_CPM_${bin_size}.bw --bigwig2 ${bam_dir}/${control}_CPM_${bin_size}.bw \
 -o "${bam_dir}/log2_CPM_${treatment}_${bin_size}.bw"

module unload bwa
module unload samtools
module unload ucsctools
module unload deeptools
date
