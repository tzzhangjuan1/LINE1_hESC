split_fastq() {
  src_dir=$1
  fastq_filename=$2
  dest_dir=$3
  #
  #cp "$src_dir/$fastq_filename" "$dest_dir"
  #paste - - - - < test.fq \
  #| tee >(awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 1:N")) print $1,$2,$3,$4}' > test.r1.fq ) \
  #| awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 2:N")) print $1,$2,$3,$4}' > test.r2.fq 
  zcat "$src_dir/${fastq_filename}.fastq.gz" | paste - - - - \
  | tee >(awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 1:N")) print $1,$2,$3,$4}' > "$dest_dir/${fastq_filename}.R1.fastq" ) \
  | awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 2:N")) print $1,$2,$3,$4}' > "$dest_dir/${fastq_filename}.R2.fastq" 
  #
  gzip "$dest_dir/${fastq_filename}.R1.fastq"
  gzip "$dest_dir/${fastq_filename}.R2.fastq"
}
filenames=(s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12)
for i in ${!filenames[@]};
do
  split_fastq data "${filenames[$i]}" NAD_data_split
done
