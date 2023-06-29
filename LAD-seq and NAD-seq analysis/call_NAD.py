import csv
import pandas as pd
import numpy as np
import pybedtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statistics import mean

def format_NAD(in_fn, out_fn):
  LAD_df = pd.read_csv(in_fn, sep='\t', skiprows=1, names=["chrom", "start", "end", "pvalue", "score", "strand", "ChIPCount", "InputCount","FDR", "log2FoldChange"])
  LAD_df[["chrom", "start", "end"]].to_csv(out_fn, index=False, header=False, sep='\t')

def call_NAD(file1, file2, file3, tmp_dir):
  out_fn_1 = tmp_dir + "/f1.tmp"
  out_fn_2 = tmp_dir + "/f2.tmp"
  out_fn_3 = tmp_dir + "/f3.tmp"

  format_NAD(file1, out_fn_1)
  format_NAD(file2, out_fn_2)
  format_NAD(file3, out_fn_3)
  int_df_12 = pybedtools.BedTool(out_fn_1)\
                       .intersect(b=out_fn_2)\
                       .to_dataframe()
  int_df_13 = pybedtools.BedTool(out_fn_1)\
                       .intersect(b=out_fn_3)\
                       .to_dataframe()
  int_df_23 = pybedtools.BedTool(out_fn_2)\
                       .intersect(b=out_fn_3)\
                       .to_dataframe()

  all_df = pd.concat([int_df_12, int_df_13, int_df_23])
  merge_df = pybedtools.BedTool.from_dataframe(all_df).sort().merge(d = 10000).to_dataframe()
  merge_df.columns=["chrom", "start", "end"]
  merge_df["length"] = merge_df.end - merge_df.start
  merge_df = merge_df[merge_df.length >= 10000][["chrom", "start", "end"]]
  return merge_df

call_NAD("results/NAD_naive_r1_10000_5.NAD.bed",\
         "results/NAD_naive_r2_10000_5.NAD.bed",\
         "results/NAD_naive_r3_10000_5.NAD.bed",\
         "results").to_csv("results/naive_NAD.bed", index=False, header=False, sep='\t')
call_NAD("results/NAD_primed_r1_10000_5.NAD.bed",\
         "results/NAD_primed_r2_10000_5.NAD.bed",\
         "results/NAD_primed_r3_10000_5.NAD.bed",\
         "results").to_csv("results/primed_NAD.bed", index=False, header=False, sep='\t')

