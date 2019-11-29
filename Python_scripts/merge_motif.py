import pandas as pd, numpy as np
from glob import glob
import re, os

#################################
# Merge with filtered motifs:
output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"

# Motif co-occurrence data:
read_df = pd.read_csv("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cooccurrence_analysis/dbf_tf/final_cobinding_motif_perc_overlap_sorted_dbf_tf.bed", sep='\t')
read_df = read_df.iloc[:,1:]
read_df["tfmotif_id"] = read_df["tfmotif_id"].str.replace(r"\|MOTIF", ".").str.replace("_human", "").str.replace("\[FLAG\]", "").str.upper()
read_df["cobind_tfmotif_id"] = read_df["cobind_tfmotif_id"].str.replace(r"\|MOTIF", ".").str.replace("_human", "").str.replace("\[FLAG\]", "").str.upper()

# Filtered motif:
filt_df = pd.read_csv(join(output_dir, "Motif_factor_and_regulatory_region_association.txt"), sep="\t")
filt_df["tf_motif_id"] = filt_df["tf_motif_id"].str.upper()
merged_df1 = pd.merge(filt_df, read_df, left_on="tf_motif_id", right_on="tfmotif_id", how="left")
merged_df = pd.merge(filt_df, read_df, left_on="tf_motif_id", right_on="tfmotif_id")

# Motifs with no motif co-occurrence:
lone_wolf_df = merged_df1[~merged_df1["tf_motif_id"].isin(merged_df["tf_motif_id"])]
lone_wolf_df["tf_motif_id"].shape

# Filter cobinding motif using main motif list:
merged_df_final = merged_df.loc[merged_df["cobind_tfmotif_id"].isin(filt_df["tf_motif_id"])]
merged_df_final_5 = merged_df_final.loc[merged_df_final["percent_overlap"] > 5]
motif_cooccur_df = merged_df_final.pivot(index="tf_motif_id", columns="cobind_tfmotif_id", values="percent_overlap")
motif_cooccur_df_5 = merged_df_final_5.pivot(index="tf_motif_id", columns="cobind_tfmotif_id", values="percent_overlap")
motif_cooccur_df.to_csv(join(output_dir,"final_motif_cooccurrence_perc_overlap_heatmap_data.txt"), sep="\t", header=True, index=True)
motif_cooccur_df_5.to_csv(join(output_dir,"final_motif_cooccurrence_perc_overlap_heatmap_data_5percentmore.txt"), sep="\t", header=True, index=True)


