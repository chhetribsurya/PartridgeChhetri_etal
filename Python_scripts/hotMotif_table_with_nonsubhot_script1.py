#!/usr/bin/env python
import pandas as pd
import os, re, collections
import pybedtools
from glob import glob
from os.path import join, splitext
from os.path import basename, dirname  

#sample = "HepG2"
sample = "HepG2full"

main_dir = "/gpfs/gpfs1/home/schhetri/for_ENCODE_jill"
main_dir = join("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill", sample) 
hot_fimomotif_file="Hot_fimomotif_intersect.filt.bed.txt"

output_dir = join(main_dir, "cisbpmotif_id_extract")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Read hotmotif file with motif intersect to compute hotmotif fraction:
# df_hotfile = pd.read_csv(join(dir_name, "MasterBedFile.HepG2.bed.srt.mrg.filt." + index_id + ".hot.final"), sep="\t", header=None)
# total_hotcount = df_hotfile.shape[0]

# Files containing HOT sites and fimo motif intersect:
# ls /gpfs/gpfs1/home/schhetri/for_ENCODE_jill/HepG2/permutation_hotsite_analysis/for_tfcount_*/permutation_num_*/Hot_fimomotif_intersect* > Hot_fimomotif_intersect_comprehensive_filelist.txt
# tf_files = open("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/HepG2_Hotsites_fimomotif_intersect_filelist.txt", "rb")

## Jill ENCODE project:
#hotmotif_filepath = "/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/cisbpmotif_id_extract/Hot_fimomotif_intersect.78.bed.txt"
hotmotif_filepath = join("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill", sample, "cisbpmotif_id_extract", hot_fimomotif_file)

dir_name = dirname(hotmotif_filepath)
base_name = basename(hotmotif_filepath)
#tf_count = base_name.split(".")[-3]
#permute_count = 1
#index_id = tf_count 
print("\nProcessing : {}".format(base_name))

# Read final hotmotif file with HOT motif intersection:
df_read = pd.read_csv(hotmotif_filepath, sep="\t", header=None)
df_read.columns = ["chrom", "start", "end", "tf_bound", "tf_count", "motif_chr", "motif_start", "motif_end", "count", "motif_bound", "motif_count"]

df_read["tf_bound"] = df_read["tf_bound"].str.replace(r"_\d+", "").str.replace("\[FLAG\]", "").str.replace(r"_human|_iso1|_iso2|_v1|_v2", "").str.upper()
df_read["motif_bound"] = df_read["motif_bound"].str.replace(r"_.*?\)", "").str.replace("\(", "").str.replace("\)", "").str.replace(r"var", "").str.replace(r"\..*", "").str.upper()  

df_motif_sum = df_read.groupby(["chrom", "start", "end", "tf_bound", "tf_count"]).apply(lambda x: x["motif_count"].sum()).reset_index(name="motif_count")
df_motif_id = df_read.groupby(["chrom", "start", "end", "tf_bound", "tf_count"]).apply(lambda x : ",".join(x["motif_bound"])).reset_index(name="motif_bound")
df_final = pd.merge(df_motif_sum, df_motif_id)

# Find uniq count of direct motifs:
df_final["motif_bound_uniqlist"] = df_final.apply(lambda x : list(set(x["motif_bound"].split(","))), axis=1)
df_final["motif_count_uniq"] = df_final.apply(lambda x : len(x["motif_bound_uniqlist"]), axis=1)

df_final["tf_bound_uniq"] = df_final["tf_bound"].str.replace(r"_\d+", "")
df_final["tf_bound_uniqlist"] = df_final["tf_bound_uniq"].apply(lambda x : list(set(x.split(","))))
df_final["tf_count_uniq"] = df_final["tf_bound_uniq"].apply(lambda x : len(set(x.split(","))))
all(df_final["tf_count"] == df_final["tf_count_uniq"])

# Filter the hot motif list table:
select_cols = ["chrom", "start", "end", "tf_bound_uniqlist", "motif_bound_uniqlist", "tf_count_uniq", "motif_count_uniq"]
df_final_filt = df_final.loc[:,select_cols]
df_final_filt["actual_motif_bound_uniqlist"] = df_final.apply(lambda x : [eachmotif for eachmotif in x["motif_bound_uniqlist"] if eachmotif in x["tf_bound_uniqlist"]], axis=1)
df_final_filt["actual_motif_count_uniq"] = df_final_filt.apply(lambda x : len(x["actual_motif_bound_uniqlist"]), axis=1)
select_cols = ["chrom", "start", "end", "tf_bound_uniqlist", "actual_motif_bound_uniqlist", "tf_count_uniq", "actual_motif_count_uniq"]
df_final_filt = df_final_filt.loc[:,select_cols]
df_final_filt.columns = ["chrom", "start", "end", "hotTF_bound_uniq", "hotmotif_present_uniq", "hotTF_bound_uniqcount", "hotmotif_present_uniqcount"]

# Hotmotif table with actual motifs or at least 1 motif:
df_motif_one = df_final_filt[df_final_filt["hotmotif_present_uniqcount"] >= 1]
df_motif_one.to_excel(join(dir_name, "Hotmotif_table_cln_min1motif.xls"), header=True, index=False)
df_motif_one.to_csv(join(dir_name, "Hotmotif_table_cln_min1motif.txt"), sep="\t", header=True, index=False)

# Table with motif ideally present; however, 0 motif presnt for corresponding TFs at HOT sites
df_motif_zero1 = df_final_filt[df_final_filt["hotmotif_present_uniqcount"] < 1]

# Read outersect hotmotif table with zero motifs to combine with absent motif table: 
hotmotif_outersect_filepath = hotmotif_filepath.rstrip().replace("intersect", "outersect")
df0_read = pd.read_csv(hotmotif_outersect_filepath, sep="\t", header=None)
df0_read.columns = ["chrom", "start", "end", "hotTF_bound_uniq", "tf_count"]    

#df0_read["tf_bound_uniq"] = df0_read["hotTF_bound_uniq"].str.replace("\[FLAG\]", "").str.replace("_human|_v.*|_iso.*", "").str.upper()  
df0_read["tf_bound_uniq"] = df0_read["hotTF_bound_uniq"].str.replace(r"_\d+", "").str.replace("\[FLAG\]", "").str.replace(r"_human|_iso1|_iso2|_v1|_v2", "").str.upper()
df0_read["tf_bound_uniqlist"] = df0_read["tf_bound_uniq"].apply(lambda x : list(set(x.split(","))))
df0_read["tf_count_uniq"] = df0_read["tf_bound_uniq"].apply(lambda x : len(set(x.split(","))))
all(df0_read["tf_count"] == df0_read["tf_count_uniq"])

# Filter the hot motif list table:
select_cols = ["chrom", "start", "end", "tf_bound_uniqlist", "motif_bound_uniqlist", "tf_count_uniq", "motif_count_uniq"]
df0_read = df0_read.loc[:,select_cols]
df0_read["motif_count_uniq"] = 0
df0_read["motif_bound_uniqlist"] = df0_read["motif_bound_uniqlist"].apply(lambda x: [])
df0_read.columns = ["chrom", "start", "end", "hotTF_bound_uniq", "hotmotif_present_uniq", "hotTF_bound_uniqcount", "hotmotif_present_uniqcount"]

concat_list = [df_motif_one, df_motif_zero1, df0_read]
combined_df = pd.concat(concat_list, ignore_index=True)
combined_df["hotTF_bound_uniqcount"] = combined_df["hotTF_bound_uniqcount"].astype(int)
combined_df["hotmotif_present_uniqcount"] = combined_df["hotmotif_present_uniqcount"].astype(int)
combined_df.to_excel(join(dir_name, "Hotmotif_table_cln_allhotsites_min0motif.xls"), header=True, index=False)
combined_df.to_csv(join(dir_name, "Hotmotif_table_cln_allhotsites_min0motif.txt"), sep="\t", header=True, index=False)

## Note : list converted to string while saving original df containing list:
# from ast import literal_eval
# final_hotmotif_df1 = pd.read_csv("/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/cisbpmotif_id_extract/Hotmotif_table_cln_allhotsites_min0motif.txt", sep="\t")
# final_hotmotif_df1["hotTF_bound_uniq"] = final_hotmotif_df1["hotTF_bound_uniq"].apply(literal_eval)
# final_hotmotif_df1["hotmotif_present_uniq"] = final_hotmotif_df1["hotmotif_present_uniq"].apply(literal_eval)


#########################################################

# 100 random permutation for HOTMOTIF Permutation plot:

#########################################################

# Final hotmotif table dataframe from previous result:
final_hotmotif_df = combined_df.copy()

# Generate uniq tf list:
total_hotcount = final_hotmotif_df.shape[0]
uniq_tf_list = []
for each_list in final_hotmotif_df["hotTF_bound_uniq"]:
    for tf in each_list:
        if tf not in uniq_tf_list:
            uniq_tf_list.append(tf)

tf_countlist = range(2,12,2) #135
tf_countlist = range(12,165,3) #135
random_seeds = range(1,101)
#tf_countlist = range(6,10,2) #135
#random_seeds = range(1,5)

tf_counts = []
permute_counts = []
true_hotsite_percent = []
one_motif_percent = []
two_motif_percent = []
three_motif_percent = []
zero_motif_percent = []

import random, time
start_time = time.time()
for tf_count in tf_countlist:
    for each_seed in random_seeds:
        #tf_count = 100
        #each_seed = 5
        tf_counts.append(tf_count)
        permute_counts.append(each_seed)
        final_hotmotif_df = combined_df.copy() # make a fresh copy of df
        random.seed(each_seed)
        sample_list = random.sample(uniq_tf_list, tf_count)
        print("Processing tf_count {} : permutation {}".format(tf_count, each_seed))
        # print(sample_list)

        final_hotmotif_df["tf_overlap"] = final_hotmotif_df.apply(lambda x: set(sample_list) & set(x["hotTF_bound_uniq"]), axis=1) 
        final_hotmotif_df["tf_overlap"] = final_hotmotif_df["tf_overlap"].apply(list)
        final_hotmotif_df["motif_overlap"] = final_hotmotif_df.apply(lambda x: set(x["tf_overlap"]) & set(x["hotmotif_present_uniq"]), axis=1) 
        final_hotmotif_df["motif_overlap"] = final_hotmotif_df["motif_overlap"].apply(list)
        
        final_hotmotif_df["tf_overlap_count"] = final_hotmotif_df.apply(lambda x : len(x["tf_overlap"]), axis=1)
        final_hotmotif_df["motif_overlap_count"] = final_hotmotif_df.apply(lambda x : len(x["motif_overlap"]), axis=1)

        # True hotsite dataframe filtering:
        #hotsite_markcount = (tf_count/2)
        hotsite_markcount = (tf_count/3)
        hotsite_df = final_hotmotif_df.loc[final_hotmotif_df["tf_overlap_count"] >= hotsite_markcount]
        total_hotsitecount = hotsite_df.shape[0]
        true_hotsite_perc = round((hotsite_df.shape[0]/float(total_hotcount))*100, 2)
        true_hotsite_percent.append(true_hotsite_perc)

        # compute percentage overlap with respect to motif counts:
        df_motif_zero = hotsite_df.loc[hotsite_df["motif_overlap_count"]==0]
        df_motif_one = hotsite_df.loc[hotsite_df["motif_overlap_count"]>=1]
        df_motif_two = hotsite_df.loc[hotsite_df["motif_overlap_count"]>=2]
        df_motif_three = hotsite_df.loc[hotsite_df["motif_overlap_count"]>=3]

        zero_motif_perc = round((df_motif_zero.shape[0]/float(total_hotsitecount))*100, 2)
        zero_motif_percent.append(zero_motif_perc)
        one_motif_perc = round((df_motif_one.shape[0]/float(total_hotsitecount))*100, 2)
        one_motif_percent.append(one_motif_perc)
        two_motif_perc = round((df_motif_two.shape[0]/float(total_hotsitecount))*100, 2)
        two_motif_percent.append(two_motif_perc)
        three_motif_perc = round((df_motif_three.shape[0]/float(total_hotsitecount))*100, 2)
        three_motif_percent.append(three_motif_perc)

df_final_perc = pd.DataFrame({"tf_count" : tf_counts,
            "permute_count" : permute_counts,
            "one_motif_perc" : one_motif_percent,
            "two_motif_perc" : two_motif_percent,
            "three_motif_perc" : three_motif_percent,
            "zero_motif_perc" : zero_motif_percent,
            "true_hotsite_perc" : true_hotsite_percent,
            })

df_final_ordered = df_final_perc.loc[:,["tf_count", "permute_count", "true_hotsite_perc", "zero_motif_perc", "one_motif_perc", 
                                        "two_motif_perc", "three_motif_perc"]]

df_final_ordered[["tf_count", "permute_count"]] = df_final_ordered.loc[:,["tf_count", "permute_count"]].astype(int)
df_final_ordered[["zero_motif_perc", "one_motif_perc", "two_motif_perc", "three_motif_perc"]] = df_final_ordered.loc[:,["zero_motif_perc", "one_motif_perc", "two_motif_perc", "three_motif_perc"]].astype(float)
df_final_ordered["true_hotsite_perc"] = df_final_ordered["true_hotsite_perc"].astype(float)
df_final_sorted = df_final_ordered.sort_values(["tf_count", "permute_count"]).reset_index(drop=True)
df_final_sorted.to_csv(join(output_dir, "Hot_permutation_motif_percent.final.redo.srt.filt.txt"), sep="\t", header=True, index=False)
print("Time for analysis = {}".format(time.time()-start_time))

