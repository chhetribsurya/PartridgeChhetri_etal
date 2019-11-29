import pandas as pd
import os, re, collections
import pybedtools
from glob import glob
from os.path import join, splitext
from os.path import basename, dirname  

sample = "HepG2full"
#sample = "HepG2"

main_dir = join("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill", sample)
tf_file = join("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill", sample, sample + "_tf_files_list.txt")

output_dir = join(main_dir, "cisbpmotif_id_extract")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#output_dir = join(main_dir, "permutation_hotsite_analysis")
#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)

# Find identifier(tf_name) for each CisBP id:
#infile = open("/gpfs/gpfs1/home/schhetri/Tools/final_custom_CisBP_database/Homosapiens_custom_cisbp_allmotifs.meme", "rb")
#infile = open("/gpfs/gpfs1/home/schhetri/for_jacob/quartile_motif_analysis/Homosapiens_combined_database_CISBP_JASPAR.meme", "rb")
infile = open("/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens_edit.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        #print(line)

cisbp_df = pd.Series(motif_list).str.split(expand=True)
cisbp_df.columns = ["motif", "cisbp_motif_id", "tf_name"]
cisbp_df["tf_name_final"] = cisbp_df["tf_name"].str.replace(r"_.*", "").str.replace("\(", "").str.replace("\)", "").str.replace(r"var", "").str.replace(r"\..*", "").str.upper()
#cisbp_df[["cisbp_motif_id", "cisbp_motif_id1"]]  = cisbp_df["cisbp_motif_id"].str.split("@").apply(pd.Series)

# Merge CisBP id with 208 HepG2 TF names:
df_read = pd.read_csv(tf_file, header=None)
df_read.columns = ["TF_name"]
df_read["tf_name_final"] = df_read["TF_name"].str.replace("\[FLAG\]", "").str.replace(r"_human|_iso1|_iso2|_v1|_v2", "").str.upper()
df_final = df_read.sort_values(["tf_name_final"]).drop_duplicates(["tf_name_final"]).reset_index(drop=True)

merged_df = pd.merge(df_final.loc[:,["tf_name_final"]], cisbp_df.loc[:,["cisbp_motif_id", "tf_name_final"]], on=["tf_name_final"], how="left")
merged_df = merged_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)
merged_df.to_csv(join(output_dir, "tf_list_with_and_withoutCisBP_id.txt"), sep="\t", header=False, index=False)
merged_df_null = merged_df.loc[merged_df["cisbp_motif_id"].isnull()]
merged_df_null.to_csv(join(output_dir, "tf_list_withoutCisBP_id.txt"), sep="\t", header=False, index=False)
merged_df.dropna(subset=["cisbp_motif_id"], inplace=True)
merged_df.to_csv(join(output_dir, "tf_list_withCisBP_id.txt"), sep="\t", header=False, index=False)

## Similarly find available motif id for all randomly permuted factors for HOT sites: 
## First generate file list with full path using bash : 
## Usage: readlink -f HepG2/permutation_hotsite_analysis/for_tfcount_*/permutation_num_*/*filelist.txt > permuted_tfs_comprehensive_filelist.txt
#string_cmd = "readlink -f " + join(main_dir, "permutation_hotsite_analysis/for_tfcount_*/permutation_num_*/*filelist.txt") + " > " + join(main_dir, "permuted_tfs_comprehensive_filelist.txt") 
#os.system(string_cmd)
#
##tf_files = open("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/permuted_tfs_comprehensive_filelist.txt", "rb")
#tf_files = open(join(main_dir, "permuted_tfs_comprehensive_filelist.txt"), "rb")
#for line in tf_files.readlines():
#    tf_filepath = line.rstrip()
#    dir_name = dirname(tf_filepath)
#    base_name = basename(tf_filepath) 
#    print("\nProcessing : {}".format(base_name))
#
#    # Merge CisBP id with 171 HepG2 TF names:
#    df_read = pd.read_csv(tf_filepath, sep="\t", header=None)
#    df_read.columns = ["TF_name"]
#    df_read["tf_name_final"] = df_read["TF_name"].str.replace(".*SL.*narrowPeak_", "").str.replace("\[FLAG\]", "").str.replace(r"_human|_iso1|_iso2|_v1|_v2", "").str.upper()
#    df_final = df_read.sort_values(["tf_name_final"]).drop_duplicates(["tf_name_final"]).reset_index(drop=True)
#
#    merged_df = pd.merge(df_final.loc[:,["tf_name_final"]], cisbp_df.loc[:,["cisbp_motif_id", "tf_name_final"]], on=["tf_name_final"])
#    merged_df = merged_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)
#    merged_df.to_csv(join(dir_name, splitext(base_name)[0] + ".motifid.txt"), sep="\t", header=False, index=False)
#
#string_cmd = "readlink -f " + join(main_dir, "permutation_hotsite_analysis/for_tfcount_*/permutation_num_*/*filelist.motifid.txt") + " > " + join(main_dir, "permuted_tfs_comprehensive_filelist.motifid.txt") 
#os.system(string_cmd)
#
### Previous line of codes deposits cisbp motif id for the given set of permuted TFs:
## Get motif pwms using memes-getmotif and scan the hotsites to find percentage of x/2 hotsites
## containing instances of >=1, >=2, >=3, motifs
## Script written in bash to do those jobs
## readlink -f HepG2/permutation_hotsite_analysis/for_tfcount_*/permutation_num_*/*filelist.motifid.txt > permuted_tfs_comprehensive_filelist.motifid.txt
#
