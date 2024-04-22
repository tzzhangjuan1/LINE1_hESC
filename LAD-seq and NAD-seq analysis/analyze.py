import csv
import pandas as pd
import numpy as np
import pybedtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
from statistics import mean
import pyBigWig
#
# Find all bins that signal > 0 of at least min_rep_num replicates
#
def call_domain(bw_filename, out_bed_filename):
  print(bw_filename)
  bw = pyBigWig.open(bw_filename)
  out_file = open(out_bed_filename, "w")
  for i in range(1,24,1):
    chrom = "chr" + str(i)
    if i == 23:
      chrom = "chrX"
    chrom_signal_vals = bw.intervals(chrom)
    for j in range(len(chrom_signal_vals)):
      if chrom_signal_vals[j][2] > 0:
         out_file.write(chrom + "\t" + str(chrom_signal_vals[j][0]) + "\t" + str(chrom_signal_vals[j][1]) + "\n")
    # end for j
  # end for i
  out_file.close()
  bw.close()  
# end call_domain
#################
def call_final_domain(fn1, fn2, fn3, out_fn):
  int_df_12 = pybedtools.BedTool(fn1)\
                       .intersect(b=fn2)\
                       .to_dataframe()
  int_df_13 = pybedtools.BedTool(fn1)\
                       .intersect(b=fn3)\
                       .to_dataframe()
  int_df_23 = pybedtools.BedTool(fn2)\
                       .intersect(b=fn3)\
                       .to_dataframe()

  all_df = pd.concat([int_df_12, int_df_13, int_df_23])
  merge_df = pybedtools.BedTool.from_dataframe(all_df).sort().merge(d = 10).to_dataframe()
  merge_df.columns=["chrom", "start", "end"]
  #merge_df["length"] = merge_df.end - merge_df.start
  #merge_df = merge_df[merge_df.length >= 10000][["chrom", "start", "end"]]
  merge_df.to_csv(out_fn, index=False, header=False, sep='\t')
# end call_final_domain
def get_domain_length(bed_fn):
  domain_df = pd.read_csv(bed_fn, sep='\t', names = ["chrom", "start", "end"])
  domain_df["length"] = domain_df["end"] - domain_df["start"]
  domain_length_dict = {"chrX": 0}
  for i in range(22):
    domain_length_dict["chr" + str(i+1)] = 0
  for i in range(len(domain_df)):
    if domain_df.loc[i, "chrom"] in domain_length_dict:
      domain_length_dict[domain_df.loc[i, "chrom"]] += domain_df.loc[i, "end"] - domain_df.loc[i, "start"]
  print(domain_length_dict)
  print(sum(domain_df["length"]))
#
def print_LAD_NAD_stat(label, LAD_fn, NAD_fn):
  print("========================= " + label)
  tmp_fn = "tmp.bed"
  pybedtools.BedTool(LAD_fn)\
            .intersect(b=NAD_fn)\
            .saveas(tmp_fn)  
  #
  print("Both LAD/NAD")
  get_domain_length(tmp_fn)
  #
  pybedtools.BedTool(LAD_fn)\
            .subtract(b=NAD_fn)\
            .saveas(tmp_fn)
  #
  print("LAD only")
  get_domain_length(tmp_fn)
  #
  pybedtools.BedTool(NAD_fn)\
            .subtract(b=LAD_fn)\
            .saveas(tmp_fn)
  #
  print("NAD only")
  get_domain_length(tmp_fn)
  #
  pybedtools.BedTool(LAD_fn)\
            .complement(genome="hg38")\
            .subtract(pybedtools.BedTool(NAD_fn))\
            .saveas(tmp_fn)
  #
  print("None of LAD/NAD")
  get_domain_length(tmp_fn)
  #
#
simple_LAD_fn = "../../results/2021_11_30_simple_NAD_LAD_call/simple_primed_LAD.bed"
epic_LAD_fn = "/mnt/work1/users/hoffmangroup/lhuynh/LINE-1/exp/2021_07_20_LAD_analysis/results/sorted_primed_LAD.bed"
simple_NAD_fn = "../../results/2021_11_30_simple_NAD_LAD_call/simple_primed_NAD.bed"
epic_NAD_fn = "/mnt/work1/users/hoffmangroup/lhuynh/LINE-1/exp/2021_08_06_NAD/results/sorted_primed_NAD.bed"
#
print("Simple primed LAD")
get_domain_length(simple_LAD_fn)
print("Epic2 primed LAD")
get_domain_length(epic_LAD_fn)
print("Simple primed NAD")
get_domain_length(simple_NAD_fn)
print("Epic2 primed NAD")
get_domain_length(epic_NAD_fn)
print_LAD_NAD_stat("SIMPLE", simple_LAD_fn, simple_NAD_fn)
print_LAD_NAD_stat("EPIC2", epic_LAD_fn, epic_NAD_fn)
#######################
#
out_dir = "../../results/2023_12_13_simple_NAD_LAD_analysis/"
naive_LAD_fn = out_dir + "simple_naive_LAD.bed"
primed_LAD_fn = out_dir + "simple_primed_LAD.bed"
naive_NAD_fn = out_dir + "simple_naive_NAD.bed"
primed_NAD_fn = out_dir + "simple_primed_NAD.bed"
primed_both_LAD_NAD_fn = out_dir + "simple_primed_LAD_NAD.bed"
primed_LAD_only_fn = out_dir + "simple_primed_LAD_only.bed"
primed_NAD_only_fn = out_dir + "simple_primed_NAD_only.bed"
primed_none_LAD_NAD_fn = out_dir + "simple_primed_none_LAD_NAD.bed"
"""
# Call NADs
NAD_bw_dir = "../../../LINE-1/exp/2021_06_01_NAD_reproducibility/results/"
for i in [1,3,5,7,9,11]:
  call_domain(NAD_bw_dir + "log2_s" + str(i) + "_10000.bw", out_dir + "s" + str(i) + ".bed")
#
call_final_domain(out_dir + "s1.bed", out_dir + "s3.bed", out_dir + "s5.bed", naive_NAD_fn) 
call_final_domain(out_dir + "s7.bed", out_dir + "s9.bed", out_dir + "s11.bed", primed_NAD_fn) 
# Call LADs
LAD_bw_dir = "../../../LINE-1/LAD_analysis/naive_pipeline/LAD_bam_files/"
bw_fns = ["s2_S2", "s5_S5", "s8_S8", "s11_S10", "s14_S13", "s17_S16"]
for fn in bw_fns:
  call_domain(LAD_bw_dir + "log2_" + fn + "_10000.bw", out_dir + fn + ".bed")
call_final_domain(out_dir + "s2_S2.bed", out_dir + "s5_S5.bed", out_dir + "s8_S8.bed", primed_LAD_fn)
call_final_domain(out_dir + "s11_S10.bed", out_dir + "s14_S13.bed", out_dir + "s17_S16.bed", naive_LAD_fn)
#
pybedtools.BedTool(primed_LAD_fn)\
          .intersect(b=primed_NAD_fn)\
          .saveas(primed_both_LAD_NAD_fn)
#
pybedtools.BedTool(primed_LAD_fn)\
          .subtract(b=primed_NAD_fn)\
          .saveas(primed_LAD_only_fn)
#
pybedtools.BedTool(primed_NAD_fn)\
          .subtract(b=primed_LAD_fn)\
          .saveas(primed_NAD_only_fn)
#
pybedtools.BedTool(primed_LAD_fn)\
          .complement(genome="hg38")\
          .subtract(pybedtools.BedTool(primed_NAD_fn))\
          .saveas(primed_none_LAD_NAD_fn)
#
"""
#get_domain_length(out_dir + "s1.bed")
#get_domain_length(out_dir + "s3.bed")
#get_domain_length(out_dir + "s5.bed")
#get_domain_length(out_dir + "s7.bed")
#get_domain_length(out_dir + "s9.bed")
#get_domain_length(out_dir + "s11.bed")
print("Naive LAD")
get_domain_length(naive_LAD_fn)
print("Primed LAD")
get_domain_length(primed_LAD_fn)
print("Naive NAD")
get_domain_length(naive_NAD_fn)
print("Primed NAD")
get_domain_length(primed_NAD_fn)
print("Both LAD/NAD (primed)")
get_domain_length(primed_both_LAD_NAD_fn)
print("LAD only( primed)")
get_domain_length(primed_LAD_only_fn)
print("NAD only (primed)")
get_domain_length(primed_NAD_only_fn)
print("None of LAD/NAD (primed)")
get_domain_length(primed_none_LAD_NAD_fn)
