import pandas as pd, numpy as np
from os.path import join, basename, splitext
from glob import glob
import subprocess
import os, shutil

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis/redo_CisBP_pwms_relabel"
#output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/novel_motif_analysis_for_3a_3b"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sub_output_dir = join(output_dir, "meme_format_motifDB")
if not os.path.exists(sub_output_dir):
    os.makedirs(sub_output_dir)

# Cisbp motif database containing PWMS:
cisbp_pwm_db_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis/CisBP_pwms"
#cisbp_pwm_db_dir = "/gpfs/gpfs1/home/schhetri/Tools/custom_CisBP_database/pwms"
cisbp_pwm_db = glob(join(cisbp_pwm_db_dir, "*.txt"))

# TF information:
cisbp_tf_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis/TF_Information_all_motifs.txt", sep="\t")
cisbp_tf_df = pd.read_csv("/gpfs/gpfs1/home/schhetri/Tools/custom_CisBP_database/TF_Information_all_motifs.txt", sep="\t")
cisbp_tf_df.columns

# Select specific cols of interest:
human_tfs = cisbp_tf_df.loc[cisbp_tf_df["TF_Species"] == "Homo_sapiens"]
human_tfs_df = human_tfs.iloc[:,range(0,13)+[16,17,18]]

# Appending the type of experminent to TF names for easy tracing later:
human_tfs_df["custom_motif_id"] = human_tfs_df["Motif_ID"] + ".txt"
human_tfs_df["custom_tf_name"] = human_tfs_df["TF_Name"] + "." + human_tfs_df["MSource_Type"]
human_tfs_df["custom_motif_id_tf_name"] = human_tfs_df["Motif_ID"] + "@" + human_tfs_df["custom_tf_name"]
human_tfs_df.head()
human_tfs_df.shape

# List comprehension of pwms:
cisbp_pwm_db_namelist = [basename(each) for each in cisbp_pwm_db]
cisbp_pwm_db_namelist[0:2]

# Filtering the motif dataframe that has motif in the database:
human_tfs_df_motif = human_tfs_df.loc[(human_tfs_df["custom_motif_id"].isin(cisbp_pwm_db_namelist))]; human_tfs_df_motif.head()
human_tfs_df_motif = human_tfs_df_motif.reset_index(drop=True)
human_tfs_df_motif = human_tfs_df_motif.sort_values(["Cutoff"], ascending=False).reset_index(drop=True).sort_values(["TF_Status"]).reset_index(drop=True)
human_tfs_df_motif["idx_custom_motif_id"] = human_tfs_df_motif["custom_motif_id"] + "@" + map(str, human_tfs_df_motif.index.tolist())
human_tfs_df_motif["custom_motif_id_tf_name"] = human_tfs_df_motif["custom_motif_id_tf_name"] + "@" + map(str, human_tfs_df_motif.index.tolist())
human_tfs_df_nomotif = human_tfs_df.loc[~(human_tfs_df["custom_motif_id"].isin(cisbp_pwm_db_namelist))]; human_tfs_df_nomotif.head()

# Listing all the databases used: (n=8) 
unique_db_used = human_tfs_df_motif["MSource_Type"].unique()

# Generate dictionary of motif file with corresponding labels + "TF_name" + "experiment"
pwm_motif_dict = human_tfs_df_motif.loc[:,["idx_custom_motif_id", "custom_motif_id_tf_name"]].set_index("idx_custom_motif_id").to_dict()

for key, value in pwm_motif_dict["custom_motif_id_tf_name"].iteritems():
    print("\nProcessing {} : {}".format(key, value))
    src_file = join(cisbp_pwm_db_dir, key.split("@")[0])
    dest_file = join(output_dir, value +  ".txt")
    shutil.copy(src_file, dest_file)
