#!/usr/bin/env python

#########################################################################
# IDEAS piechart analysis-discrepancy impact-ratio and heatmap clustering
#########################################################################

import pandas as pd, numpy as np
from glob import glob
import re
from os.path import join, basename, dirname

file_list = glob("/gpfs/gpfs1/home/schhetri/paper_revision_3/combined_piecharts_quartile_motifs_plus_original208//*intersectbed_counts_with_ideas.txt")
#peak_file_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*")
dir_path = "/gpfs/gpfs1/home/schhetri/paper_revision_3/combined_quartile_plus_original208_bed/*srt*"
peak_file_list = glob(dir_path)

output_dir = "~/paper_revision_3"
df_list = []
tf_list = []

regex_pat = re.compile(r"_intersectbed_counts_with_ideas.txt")
for each_file in file_list:
    tf_name = regex_pat.split(basename(each_file))[0]
    df = pd.read_csv(each_file, sep="\t", header=None)
    df.columns = ["state", "peak_count"]
    df["tf_name"] = tf_name
    df = df.loc[:,["tf_name", "state", "peak_count"]]
    df_list.append(df)
    tf_list.append(tf_name)

combined_df = pd.concat(df_list, ignore_index=True)
final_df = combined_df.sort_values(["tf_name"]).reset_index(drop=True)
#sum_df = final_df.groupby(["tf_name"])["peak_count"].sum().reset_index(drop=False)
#final_df["total_count"] = final_df.groupby(["tf_name"])["peak_count"].transform(sum)
#final_df["fraction_bind"] = final_df["peak_count"].astype(float)/final_df["total_count"]
# count_df = final_df.groupby(["tf_name"])["total_count"].value_counts()

peak_tf_list = []
peak_count_list = []
#regex_pat = re.compile(r"narrowPeak_")
regex = re.compile(r".*/(.*).srt.(.*).bed")
narrow_regex = re.compile(r".*narrowPeak_(.*)$")

for each_file in peak_file_list:
    if "narrowPeak" in each_file:
        tf_name = narrow_regex.split(basename(each_file))[1]
        df = pd.read_csv(each_file, sep="\t", header=None)
        peak_count_list.append(df.shape[0])
        peak_tf_list.append(tf_name)
    else:
        capture_list = regex.findall(each_file)[0]
        tf_name = capture_list[0] + "." + capture_list[1]
        df = pd.read_csv(each_file, sep="\t", header=None)
        peak_count_list.append(df.shape[0])
        peak_tf_list.append(tf_name)

peak_df = pd.DataFrame({"tf_name": peak_tf_list, 
                        "total_peak_counts" : peak_count_list })
peak_df = peak_df.loc[:,["tf_name","total_peak_counts"]]
peak_df = peak_df.sort_values(["tf_name"]).reset_index(drop=True)
peak_df.to_csv(join(output_dir, "quartilebased_plus208_all_tfs_peakcount.txt"), sep ="\t", header = True, index = False)

merged_df = pd.merge(final_df, peak_df, on="tf_name")
merged_df["fraction_bind"] = merged_df["peak_count"]/merged_df["total_peak_counts"].astype(float)
merged_df.to_csv(join(output_dir, "quartilebased_plus208_all_tfs_merged_peakcount.txt"), sep ="\t", header = True, index = False)

heatmap_df = merged_df.pivot(index="tf_name", columns="state", values="fraction_bind")
#heatmap_df["REST_STATE"] = heatmap_df.iloc[:,1:].apply(lambda x : 1-sum(x), axis =1)
heatmap_df["Total_Frac"] = heatmap_df.iloc[:,1:].apply(lambda x : sum(x), axis =1)
heatmap_df.to_csv(join(output_dir, "quartilebased_plus208_ideastyle_promoter_state_dist_for_heatmap.txt"), sep ="\t", header = True, index = True)

#df = pd.concat([peak_df, sum_df], axis=1)
#df["diff"] = df["peak_counts"] - df["peak_count"]
#df["impact_ratio"] = df["diff"].astype(float)/df["peak_counts"].astype(float)
#
## Explaination note: All the peaks with less than 50% overlap is disqualified for peak distribution
#diff_df = df.iloc[df["diff"].nlargest(10).index]
#
## The impact ratio of disqualified peaks shows that max impact would be of CBX1, with shift of just 0.4%(0.004831*100) 
## at max, that is percentage shift it might have caused - had those disqualified peaks been included in the piechart. 
#impact_df = df.iloc[df["impact_ratio"].nlargest(10).index]
