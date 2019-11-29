
####################################################

# Motif and Factors relation with Prom/Enh/Prom-Enh

####################################################
# Run Locally:
import pandas as pd
from os.path import os, join

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
input_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"

# R :library("readxl") 
#ideas_table <- read_excel(file.path(output_dir, "IDEAS_reannotation_table.xlsx" ))
ideas_table = pd.read_excel(join(input_dir, "IDEAS_reannotation_table.xlsx"))
tomtom_motif_df = pd.read_csv(join(input_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_redo.txt"), sep="\t")

# Subset df:
tomtom_motif_df["anno_new"] = tomtom_motif_df["anno"].replace({"match":"concordance", "mismatch":"discordance"})
tomtom_subset_df = tomtom_motif_df.loc[:, ["TF_NAME", "tf_motif_id", "anno_new"]]
merged_df = pd.merge(tomtom_subset_df,ideas_table, left_on="TF_NAME", right_on="Target", indicator=True)
merged_df["final_ano"] = merged_df["annotation"].str.replace(r"Cluster .* \(", "").str.replace(r"regions\)", "").str.replace(r"\)", "")
merged_df.to_csv(join(output_dir, "Motif_factor_and_regulatory_region_association_redo.txt"), sep="\t", header=True, index=False)
