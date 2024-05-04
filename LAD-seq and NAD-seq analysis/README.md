# LAD-seq and NAD-seq analysis
## For NADs:
  Raw fastq files -> split.sh -> R1/R2 fastq files -> NAD.slurm.sh -> bam files -> visual.slurm.sh -> bw files -> analyze.py -> NADs (bed files) 
## For LADs:
  Raw fastq files -> LAD.slurm.sh -> bw files -> analyze.py -> LADs (bed files)
# Bed files called are in the Results/

Data processing and bed calling by Dr.Linh Huynh (vietlinh.huynh@gmail.com); Downstream analysis on the bed files by Dr.Juan Zhang (tzzhangjuan@gmail.com)
