#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name chipQC2

export PATH=/home/wranalab/lcaldwel/.local/bin/cutadapt:$PATH


cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_S1_L006_R1_001.fastq.gz -p tr_S1_L006_R2_001.fastq.gz S1_L006_R1_001.fastq.gz S1_L006_R2_001.fastq.gz

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_s2_L006_R1_001.fastq.gz -p tr_s2_L006_R2_001.fastq.gz s2_L006_R1_001.fastq.gz s2_L006_R2_001.fastq.gz

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_s3_L006_R1_001.fastq.gz -p tr_s3_L006_R2_001.fastq.gz s3_L006_R1_001.fastq.gz s3_L006_R2_001.fastq.gz

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_s4_L006_R1_001.fastq.gz -p tr_s4_L006_R2_001.fastq.gz s4_L006_R1_001.fastq.gz s4_L006_R2_001.fastq.gz

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o tr_s5_L006_R1_001.fastq.gz -p tr_s5_L006_R2_001.fastq.gz s5_L006_R1_001.fastq.gz s5_L006_R2_001.fastq.gz

export PATH=/home/wranalab/lcaldwel/programs/fastqc:$PATH

fastqc tr_*.gz

multiqc . --ignore "FastQC/"
multiqc . --ignore "multiqc_firstpass/"

echo "completed cutadapt and fastqc second pass with exit code $? at: `date`"
