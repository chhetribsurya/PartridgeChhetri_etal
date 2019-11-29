
import re, os, pickle
from os import makedirs
from os.path import exists, join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np

tf_svm_scoredict={"a":0.135, "b": -0.12, "c":-1.34,"d":1.00,"e":1.2,"f":0.15}
test_new={6:"a,b,c", 5:"c,d,e,f",7:"e,f,g", 8:"e,f,g", 9:"e,f,g"}

test_df_new = pd.DataFrame(test_new.items())
test_df_new.columns = ["num_tfbound", "tfname"]; test_df_new

# For hotsite problem 2:
master_list = []

# For hotsite problem 1:
num_tfbound_list = []
tfscore_list = []

# Hotsites dataframe - for problem 1 and 2 both:
for idx,row in test_df_new.iterrows():
    tfid_splitted=row["tfname"].split(",")
    num_tfbound =row["num_tfbound"]
    for each_id in tfid_splitted:
        print num_tfbound, each_id
        tfscore_list.append(tf_svm_scoredict.get(each_id))
        num_tfbound_list.append(num_tfbound)

        # For hotsite problem 2:
        master_list.append([idx,num_tfbound, tf_svm_scoredict.get(each_id), each_id])

# For hotsite problem 1:
df = pd.concat([pd.Series(num_tfbound_list), pd.Series(tfscore_list)], axis=1)
df.columns = ["tf_cobound", "svm_score"]

# For hotsite problem 2 (master list df)
master_df = pd.DataFrame(master_list)
master_df.columns = ["hotsite_idx", "tf_bound", "score", "tf_id"]
master_df.groupby(["hotsite_idx", "tf_bound"])["score"].mean()

# Handles the cases with hotsites containing more than 1 peaks or motifs 
# for same factor like FOXA3 and KDM2A at hotsite 4 and 6 in hepg2 hotsites:
df_grouped = master_df.groupby(["hotsite_idx", "tf_bound"])["score"].mean().reset_index(name="svmscore")
df_final = df_grouped.pivot(index="hotsite_idx", columns="tf_bound", values="svmscore")

# For hotsite problem 2
# Rank hotsites with TFs_svmscore or classifier value or less : ()
df_rank = df_final.rank(ascending=False,method="dense",pct=True)

# Filter all rows or hotsites with no 5% TF classifier value:
df_rank["5perc_present"] = df_rank.apply(lambda row : row > 0.7, axis=0).astype(int).apply(np.sum, axis=1)
df_rank["tf_cobound"] = df_grouped["tf_bound"]

# Dframe containing records for hotsites with 5% classifier value or less:
final_df = df_rank.loc[df_rank["5perc_present"] > 0 ]

# Total frequency of hotsites containing top 5% classifier value for any TF present:
# Frequency of hotsite or total with 5% classifier value present
hist(df_rank["tf_cobound"])

#####################################################
# df_rank.apply(lambda row : row > 0.2, axis=1)
# df_final.rank(ascending=False,method="dense",pct=True)

# Total frequency of hotsites containing top 5% classifier value for any TF present:
# df_rank.apply(lambda row : np.where(row > 0.2, "Y","N"))
df_rank.apply(lambda row : row > 0.7, axis=1).astype(int)
df_rank["5perc_present"] = df_rank.apply(lambda row : row > 0.7, axis=0).astype(int).apply(np.sum, axis=1)

# Give score of 1 if present that is (more than 1 TF with classifier value 0.05); 
# else assign 0 if no TF classifier value with 0.05 or 5% is present in hotsite index: 
# Then grouping by TF bound gives us the frequency of hotsites with TF classifier value 
# containing 0.05 or 5% classifier value for histogram plot:
df_rank["5perc_present"] = np.where(df_rank["5perc_present"] > 0, 1, 0 )
df_rank["tf_cobound"] = df_grouped["tf_bound"]

# Frequency of hotsite or total with 5% classifier value present
df_rank.groupby(["tf_cobound"])["5perc_present"].sum().reset_index(name="5perc_present")

