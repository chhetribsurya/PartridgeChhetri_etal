import pandas as pd
import os
import re
from pybedtools import BedTool
from glob import glob
from os.path import join
from os.path import splitext
from os.path import basename

input_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/gerp_analysis_50bp_up_and_down/split_chrom_intersect_hotspot_gerp_scores")
output_dir = join(input_dir, "analysis_output")

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


chromwise_gerp_analysis_files = glob(join(input_dir, "*intersect.bed"))
chromwise_gerp_analysis_files_df = []

for each in chromwise_gerp_analysis_files:
	df = pd.read_csv(each, sep="\t", header=None)
	df_groupby = df.groupby([0,1,2,3,4])[10].apply(lambda x: x.mean())
	df_groupby = df_groupby.reset_index()
	df_groupby.columns = ["chrom", "start", "end", "anno_state", "tf_count", "mean_rs_score"]
	
	bins = [1,5,10,20,30,40,50,70,100,500]
	names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
	df_groupby['binned_tf_count'] = pd.cut(df_groupby["tf_count"], bins, right=False, labels=names)
	df_groupby_sorted = df_groupby.sort_values(["tf_count"])
	chromwise_gerp_analysis_files_df.append(df_groupby_sorted)

combined_gerp_df = pd.concat(chromwise_gerp_analysis_files_df, ignore_index=True)
combined_gerp_df.to_csv(join(output_dir, "combined_chromwise_gerp_analysis_boxplot_data.txt"), sep="\t", header=True, index=False)


### Gerp elements analysis:
output_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/split_chrom_intersect_hotspot_gerp_elements")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

### Analysis for the overlap fraction per kilobases [(overlap/total_length)*100]
gerp_element_file = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/split_chrom_intersect_hotspot_gerp_elements/combined_hotspot_gerp_elements_intersect_wao_parameter.bed")
gerp_elem_df = pd.read_csv(gerp_element_file, sep="\t", header=None)
gerp_elem_df.rename(columns={13:"overlap_bases"}, inplace=True)
gerp_elem_df.rename(columns={10:"gerp_elem_sum_score"}, inplace=True)
gerp_elem_df.rename(columns={9:"gerp_elem_length"}, inplace=True)
gerp_elem_df.loc[gerp_elem_df["gerp_elem_length"] == ".", "gerp_elem_length"] = 0
gerp_elem_df.loc[gerp_elem_df["gerp_elem_sum_score"] == ".", "gerp_elem_sum_score"] = 0

gerp_elem_df["gerp_elem_score_perbase"] = gerp_elem_df["gerp_elem_sum_score"].astype(float)/gerp_elem_df["gerp_elem_length"].astype(float)
gerp_elem_df["gerp_elem_score_perbase"].fillna(0,inplace=True)
gerp_elem_df["overlap_score"] = gerp_elem_df["overlap_bases"] * gerp_elem_df["gerp_elem_score_perbase"]
#mask_gt = gerp_elem_df["overlap_bases"] > 0
#gerp_elem_df.loc[mask_gt, "overlap_bases"] = gerp_elem_df["overlap_bases"] # to adjust the bedtool 0 based coordinates: 

gerp_elem_df_grouped = gerp_elem_df.groupby([0,1,2,3,4])["overlap_bases", "overlap_score"].apply(lambda x: x.sum())
gerp_elem_df_overlap = gerp_elem_df_grouped.reset_index()
gerp_elem_df_overlap["total_length"] = gerp_elem_df_overlap.iloc[:,2].astype(int) - gerp_elem_df_overlap.iloc[:,1].astype(int)

gerp_elem_df_overlap["fraction_overlap_per_100bp"] = (gerp_elem_df_overlap["overlap_bases"]/gerp_elem_df_overlap["total_length"])*100
gerp_elem_df_overlap["fraction_overlap_per_kb"] = (gerp_elem_df_overlap["overlap_bases"]/gerp_elem_df_overlap["total_length"])*1000
gerp_elem_df_overlap["fraction_overlap_score_per_100bp"] = (gerp_elem_df_overlap["overlap_score"]/gerp_elem_df_overlap["total_length"])*100
gerp_elem_df_overlap["fraction_overlap_score_per_kb"] = (gerp_elem_df_overlap["overlap_score"]/gerp_elem_df_overlap["total_length"])*1000

bins = [1,5,10,20,30,40,50,70,100,150,500]
names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100-149", "150+" ]
gerp_elem_df_overlap['binned_tf_count_zoomed'] = pd.cut(gerp_elem_df_overlap.iloc[:,4], bins, right=False, labels=names)

bins = [1,5,10,20,30,40,50,70,100,500]
names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
gerp_elem_df_overlap['binned_tf_count'] = pd.cut(gerp_elem_df_overlap.iloc[:,4], bins, right=False, labels=names)

gerp_elem_df_sorted = gerp_elem_df_overlap.sort_values([4]).reset_index(drop=True)
gerp_elem_df_sorted.to_csv(join(output_dir, "combined_gerp_element_analysis_boxplot_data_with_gerp_elem_overlap_fraction_final.txt"), sep="\t", header=True, index=False)

