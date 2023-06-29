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
#from enrichment import analyze_enrichment
#
#
#
def analyze_permutation(domain_fn, RE_df, total_genome_length):
  RE_shuffle_df = pybedtools.BedTool\
                            .from_dataframe(RE_df)\
                            .shuffle(genome='hg38', chrom=True)\
                            .to_dataframe()
  #
  RE_shuffle_in_domain_df = pybedtools.BedTool\
                                      .from_dataframe(RE_shuffle_df)\
                                      .sort()\
                                      .intersect(b = domain_fn, wa=True)\
                                      .to_dataframe()
  #
  RE_in_domain_df = pybedtools.BedTool\
                              .from_dataframe(RE_df)\
                              .sort()\
                              .intersect(b = domain_fn, wa=True)\
                              .to_dataframe()
  return len(RE_in_domain_df), len(RE_shuffle_in_domain_df)
#
def get_domain_length(domain_fn):
  domain_df = pd.read_csv(domain_fn, sep='\t', names = ["chrom", "start", "end"])
  domain_df["length"] = domain_df["end"] - domain_df["start"]
  return domain_df["length"].sum()
  
#
def analyze_enrichment(LAD_fn, NAD_fn, L1_filter, exp_name):
  print(LAD_fn)
  print(NAD_fn)
  #
  print("Start analyzing enrichment ... ")
  print(get_domain_length(LAD_fn))
  print(get_domain_length(NAD_fn))
  out_dir = "../../results/2023_05_04_Enrichment/"
  total_genome_length = 3031042417
  #
  # Read all UCSC repetitive elements
  all_re_df = pd.read_csv("/mnt/work1/users/hoffmangroup/lhuynh/RepeatHiC/data/2020_04_07_L1_hg38_ref/All_REs.hg38.ucsc", sep='\t')
  all_re_df = all_re_df[["genoName", "genoStart", "genoEnd", "repName", "repClass", "repFamily"]]
  #
  #re_count_df = all_re_df.repName.value_counts().to_frame()
  #re_count_df.index.name = "RE_name"
  #re_count_df.reset_index(inplace=True)
  #re_count_df.columns = ["RE_name", "freq"]
  #rep_name_df = re_count_df[re_count_df.freq >= 1000]
  #
  repf_count_df = all_re_df.repFamily.value_counts().to_frame()
  repf_count_df.index.name = "RE_family"
  repf_count_df.reset_index(inplace=True)
  repf_count_df.columns = ["RE_family", "freq"]
  rep_family_df = repf_count_df[repf_count_df.freq >= 1000]
  #print(rep_family_df)
  #
  #2	Simple_repeat   714283
  #3              MIR   610603
  #4               L2   481336
  #5        ERVL-MaLR   363032
  #6      hAT-Charlie   268067
  #7             ERV1   184079
  #8             ERVL   169746
  rep_family_df = rep_family_df[rep_family_df["RE_family"].isin(["L1", "Alu", "Simple_repeat", "MIR", "L2", "ERVL-MalR", "ERV1", "ERVL"])]
  #
  rep_family_df["LAD_coverage_shuffle"] = \
  rep_family_df["NAD_coverage_shuffle"] = \
  rep_family_df["LAD_coverage"] = \
  rep_family_df["NAD_coverage"] = rep_family_df["freq"]
  #
  chrom_list = ["chr" + str(i) for i in range(1, 23)]
  chrom_list.append("chrX")
  #
  for index, row in rep_family_df.iterrows():
    RE_by_family_df = all_re_df[(all_re_df["repFamily"] == row["RE_family"]) & \
                                (all_re_df["genoName"].isin(chrom_list))][["genoName", "genoStart", "genoEnd"]]
    if row["RE_family"] == "L1":
      if L1_filter == "short":
        RE_by_family_df = all_re_df[(all_re_df["repFamily"] == row["RE_family"]) & \
                                    (all_re_df["genoName"].isin(chrom_list)) & \
                                    (all_re_df["genoEnd"] - all_re_df["genoStart"] >= 5000)][["genoName", "genoStart", "genoEnd"]]
      if L1_filter == "young":
        RE_by_family_df = all_re_df[(all_re_df["repFamily"] == row["RE_family"]) & \
                                    (all_re_df["repName"].isin(["L1PA" + str(i) for i in range(2,9)])) & \
                                    (all_re_df["genoName"].isin(chrom_list))][["genoName", "genoStart", "genoEnd"]]
    print([row["RE_family"], len(RE_by_family_df.index)])
    #print(young_L1_df[young_L1_df["RE_family"] == "L1"])
    #print([len(rep_family_df.index), len(short_L1_df.index), len(young_L1_df.index)])

    LAD_coverage, LAD_coverage_shuffle = analyze_permutation(LAD_fn, RE_by_family_df, total_genome_length)
    NAD_coverage, NAD_coverage_shuffle = analyze_permutation(NAD_fn, RE_by_family_df, total_genome_length)
    #
    rep_family_df.at[index, ["LAD_coverage", "NAD_coverage", "LAD_coverage_shuffle", "NAD_coverage_shuffle"]] \
      = (LAD_coverage, NAD_coverage, LAD_coverage_shuffle, NAD_coverage_shuffle)
  #
  rep_family_df.to_csv(out_dir + "RE_enrichment_by_family" + "_" + exp_name + ".txt", index=False, sep='\t')

  #import csv
  #import pandas as pd
  #import numpy as np
  #import matplotlib
  #matplotlib.use('Agg')
  #import matplotlib.pyplot as plt
  #import seaborn as sns
  #from scipy import stats
  #from scipy.stats import ttest_ind
  #from statistics import mean
  #
  #
  #
  #RE_family	freq	LAD_coverage_shuffle	NAD_coverage_shuffle	LAD_coverage	NAD_coverage
  #Alu	1262425	465939	554310	291769	626460
  #L1	1017024	423692	423641	450697	440607
  #enrichment_df = pd.read_csv("../../results/2023_05_04_Enrichment/RE_enrichment_by_family.txt", sep='\t')
  rep_family_df["LAD_enrichment"] = rep_family_df["LAD_coverage"]/rep_family_df["LAD_coverage_shuffle"]
  rep_family_df["NAD_enrichment"] = rep_family_df["NAD_coverage"]/rep_family_df["NAD_coverage_shuffle"]
  #
  plt.figure(figsize=(5,7))
  bar_width = 0.35
  bar_positions = np.arange(len(rep_family_df))
  # create the bar chart
  plt.bar(bar_positions-bar_width/2,rep_family_df["LAD_enrichment"], bar_width, label="LAD")
  plt.bar(bar_positions+bar_width/2,rep_family_df["NAD_enrichment"], bar_width, label="NAD")
  plt.axhline(y=1, color='black', linestyle='--')
  #
  plt.xticks(bar_positions, rep_family_df["RE_family"], rotation=45)
  plt.ylabel("Enrichment fold compared to random regions")
  plt.legend()
  #
  plt.savefig("/mnt/work1/users/hoffmangroup/www/people/lhuynh/internal/2023/LAD_NAD_RE_enrichment" + "_" + exp_name + ".png")
  plt.close()
  #
#
simple_LAD_fn = "../../results/2021_11_30_simple_NAD_LAD_call/simple_primed_LAD.bed"
epic_LAD_fn = "/mnt/work1/users/hoffmangroup/lhuynh/LINE-1/exp/2021_07_20_LAD_analysis/results/sorted_primed_LAD.bed"
ref_LAD_fn = "hESC_LAD.bed.txt"
simple_NAD_fn = "../../results/2021_11_30_simple_NAD_LAD_call/simple_primed_NAD.bed"
epic_NAD_fn = "/mnt/work1/users/hoffmangroup/lhuynh/LINE-1/exp/2021_08_06_NAD/results/sorted_primed_NAD.bed"
#
# 
analyze_enrichment(simple_LAD_fn, simple_NAD_fn, "all", "simple_caller")
analyze_enrichment(epic_LAD_fn, epic_NAD_fn, "all", "epic_caller")
analyze_enrichment(ref_LAD_fn, epic_NAD_fn, "all", "ref_LAD")
#
analyze_enrichment(epic_LAD_fn, epic_NAD_fn, "young", "young_epic_caller")
analyze_enrichment(epic_LAD_fn, epic_NAD_fn, "short", "long_epic_caller")

