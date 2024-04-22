#!/bin/bash


echo "calling peaks round 1 with macs..."
# -t is treatment.
# -n is name of outfile.

# for file in *.bam; do
macs3 callpeak --nomodel -B -f BAMPE -t filtered_S1.bam -n ctr-input --nolambda --keep-dup all --gsize hs
macs3 callpeak --nomodel -B -f BAMPE -t filtered_s2.bam -n MALAT1-CHIRP-rep1 --nolambda --keep-dup all --gsize hs
macs3 callpeak --nomodel -B -f BAMPE -t filtered_s3.bam -n L1-CHIRP-rep1 --nolambda --keep-dup all --gsize hs
macs3 callpeak --nomodel -B -f BAMPE -t filtered_s4.bam -n MALAT1-CHIRP-rep2 --nolambda --keep-dup all --gsize hs
macs3 callpeak --nomodel -B -f BAMPE -t filtered_s5.bam -n L1-CHIRP-rep2 --nolambda --keep-dup all --gsize hs

#now do merged calling
echo "merged calling, v1..."
macs3 callpeak --nomodel -B -f BAMPE -t malat1-chirp_merge.bam -n malat1-chirp_merge --nolambda --keep-dup all --gsize hs
macs3 callpeak --nomodel -B -f BAMPE -t l1-chirp_merge.bam -n l1-chirp_merge --nolambda --keep-dup all --gsize hs

# for broad cutoff: here you use the control samples as a background.
echo "broad calling..."
macs3 callpeak -t filtered_s2.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_MALAT1-CHIRP-rep1 -f BAMPE
macs3 callpeak -t filtered_s3.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_L1-CHIRP-rep1 -f BAMPE
macs3 callpeak -t filtered_s4.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_MALAT1-CHIRP-rep2 -f BAMPE
macs3 callpeak -t filtered_s5.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_L1-CHIRP-rep2 -f BAMPE


#broad merge now.
echo "broad merged calling..."
macs3 callpeak -t malat1-chirp_merge.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_malat1-chirp_merge  -f BAMPE
macs3 callpeak -t l1-chirp_merge.bam -c filtered_S1.bam --broad -g hs --broad-cutoff 0.1 -n broad_l1-chirp_merge -f BAMPE

echo "Finished chipseq MACS3 run with exit code $? at: `date`"
