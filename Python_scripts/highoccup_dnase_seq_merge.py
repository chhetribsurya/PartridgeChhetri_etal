import pandas as pd
import pybedtools
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
#from rpy2.robjects.packages import importr  #stats = importr('stats')
#from rpy2.robjects.vectors import FloatVector
from os.path import join
from os.path import basename
from os.path import splitext
import warnings


ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_dnase/merged_d_100"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

""" Dnase window size for TF intersection 50bp up and downstream from center"""
size_window = 0 # either 0 or greater than 0

""" Select the state number, to be analysed on """
#select_state_num = range(1,8+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
prom_assoc_state_num = range(1,8+1)# Add +1; range() is exclusive:
strong_enh_state_num = [15,16] # Add +1; range() is exclusive:
weak_enh_state_num = range(17,20+1) # Add +1; range() is exclusive:


""" Merging dnase 2 reps df"""
dnase_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_dnase/hepg2_dnase"
dnase_rep1_bed_file = join(dnase_file_dir, "ENCFF505SRS_rep1_dnase_seq.bed")
dnase_rep2_bed_file = join(dnase_file_dir, "ENCFF268DTI_rep2_dnase_seq.bed")

dnase_rep1_df = pd.read_csv(dnase_rep1_bed_file, sep = "\t", header = None)
dnase_rep1_df = dnase_rep1_df.iloc[:, [0, 1, 2, 9]]
dnase_rep1_df.columns = ["chrom", "start", "end", "rep1_summit"]

dnase_rep2_df = pd.read_csv(dnase_rep2_bed_file, sep = "\t", header = None)
dnase_rep2_df = dnase_rep2_df.iloc[:, [0, 1, 2, 9]]
dnase_rep2_df.columns = ["chrom", "start", "end", "rep2_summit"]

all_dfs = [dnase_rep1_df, dnase_rep2_df]
merged_dnase_df = reduce(lambda left, right: pd.merge(left, right, on=["chrom","start","end"]), all_dfs)

### Intersect Hepg2 bed file and Filter the rows with max base overlap size for duplicated rows:
hepg2_ideas_df = pd.read_csv(ideas_hepg2_file, sep="\t")
hepg2_ideas_df = hepg2_ideas_df.iloc[:, 1:5]
hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "chromStart"])

hepg2_bed_file = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
dnase_bed_file = pybedtools.BedTool.from_dataframe(merged_dnase_df)
#dnase_merged_bed_file = dnase_bed_file.merge(d=100, c=[4,5], o=["distinct","distinct"])

dnase_ideas_intersect = dnase_bed_file.intersect(hepg2_bed_file, wao=True)

### remove the feature with 0 overlap size(last column)
dnase_ideas_intersect_df = pd.read_csv(dnase_ideas_intersect.fn, sep="\t", header=None)
intersect_df =  dnase_ideas_intersect_df[dnase_ideas_intersect_df.iloc[:, -1] > 0] 

### Filter the rows with max base overlap size for duplicated rows:(column 9 here represents the size for base overlap)
duplicated_df = intersect_df[intersect_df.duplicated([0,1,2], keep=False)]
duplicated_filtered_df = duplicated_df.loc[duplicated_df.groupby([0,1,2])[9].idxmax()]
non_duplicated_df = intersect_df[~intersect_df.duplicated([0,1,2], keep=False)] # non-duplicated df
if (intersect_df.shape[0] == duplicated_df.shape[0] + non_duplicated_df.shape[0]): # 47743 = 34578 + 13165 
    print "Duplicated Rows filtered after groupby and merged successfully"

dnase_uniq_df = pd.concat([duplicated_filtered_df, non_duplicated_df], ignore_index=True)
dnase_final_df = dnase_uniq_df.iloc[:,range(0,3)+[8]]
dnase_final_df.columns = ["chrom", "start", "end", "state"]
dnase_final_df = dnase_final_df.sort_values(["chrom", "start"])
dnase_final_df = dnase_final_df.reset_index(drop=True) # reordering of the index
dnase_ideas_distribution = dnase_final_df["state"].value_counts()

dnase_ideas_distribution.to_csv(join(output_dir, "final_dnase_wholegenome_piechart_distribution_data.txt"), header=True, index=True, sep="\t")
dnase_final_df.to_csv(join(output_dir, "final_dnase_annotation_with_ideas_states.bed"), header=True, index=True, sep="\t")


""" Read the mnemonics file, and create a dictionary of states with state number ID"""
read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
state_list = read_file["Mnemonics"]
state_num_list = [ i for i in range(1,len(state_list)+1) ]
ideas_state_dict = dict(zip(state_num_list,state_list))

#target_state = [ideas_state_dict[each] for each in select_state_num] #Depends on the select_state_num variable input
prom_assoc_state = [ideas_state_dict[each] for each in prom_assoc_state_num] #Depends on the select_state_num variable input
strong_enh_state = [ideas_state_dict[each] for each in strong_enh_state_num] #Depends on the select_state_num variable input
weak_enh_state = [ideas_state_dict[each] for each in weak_enh_state_num] #Depends on the select_state_num variable input

def merge_ideas_pybed(Ideas_state_dict, Dnase_ideas_df, Select_state_num, Merged_state_anno):
    """ Depends on the select_state_num variable input"""
    target_state = [Ideas_state_dict[each] for each in Select_state_num]
    Dnase_ideas_df["state_anno"] = Merged_state_anno
    Dnase_ideas_df["strand"] = "."  
    dnase_ideas_select_df = Dnase_ideas_df[Dnase_ideas_df["state"].isin(target_state)]  

    dnase_ideas_sorted_df = dnase_ideas_select_df.sort_values(["chrom", "start"])
    dnase_ideas_sorted_df.to_csv(join(output_dir, Merged_state_anno + "_dnase_ideas_states.bed"), header=True, sep="\t")

    """ Generate pybed object"""
    dnase_ideas_pybed = pybedtools.BedTool.from_dataframe(dnase_ideas_sorted_df) 
    merge_ideas_pybed = dnase_ideas_pybed.merge(d=100, c=[4,5], o=["distinct", "distinct"])
    return(merge_ideas_pybed)

""" Merging the combined_ideas and hepg2 file by calling a function"""
# hepg2_ideas_pybed =  merge_ideas_pybed(ideas_mnemonics_file, ideas_hepg2_file, select_state_num) 
hepg2_ideas_pybed_prom =  merge_ideas_pybed(ideas_state_dict, dnase_final_df, prom_assoc_state_num, "prom_assoc") 
hepg2_ideas_pybed_sEnh =  merge_ideas_pybed(ideas_state_dict, dnase_final_df, strong_enh_state_num, "sEnh_assoc") 
hepg2_ideas_pybed_wEnh =  merge_ideas_pybed(ideas_state_dict, dnase_final_df, weak_enh_state_num, "wEnh_assoc") 

print hepg2_ideas_pybed_prom.head()
print hepg2_ideas_pybed_sEnh.head()
print hepg2_ideas_pybed_wEnh.head()

concat_list_df = []
df_line_count = []
hepg2_ideas_list = [hepg2_ideas_pybed_prom, hepg2_ideas_pybed_sEnh, hepg2_ideas_pybed_wEnh]
for each in hepg2_ideas_list:
    pybed_df = pd.read_csv(each.fn, sep="\t", header=None)
    pybed_df.columns = ["chrom", "start","end", "state", "state_anno"]
    concat_list_df.append(pybed_df)
    df_line_count.append(pybed_df.shape[0])

dnase_ideas_df = pd.concat(concat_list_df, ignore_index=True)
dnase_ideas_df = dnase_ideas_df.sort_values(["chrom", "start"])

### For 300bp up and downstream of dnase sites:
#size_window = 300
if size_window > 0:
    output_dir = join(output_dir, "dnase_window_size_" + size_window)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    dnase_ideas_df_mod = dnase_ideas_df.copy()
    dnase_ideas_df_mod["midpoint"] = (dnase_ideas_df["start"] + dnase_ideas_df["end"])/2
    dnase_ideas_df_mod["midpoint"] = dnase_ideas_df_mod["midpoint"].astype(int)
    dnase_ideas_df_mod["chromstart"] = dnase_ideas_df_mod["midpoint"] - (size_window)
    dnase_ideas_df_mod["chromend"] = dnase_ideas_df_mod["midpoint"] + (size_window)
    dnase_ideas_pybed = pybedtools.BedTool.from_dataframe(dnase_ideas_df_mod)
else:
    dnase_ideas_pybed = pybedtools.BedTool.from_dataframe(dnase_ideas_df)


if sum(df_line_count) == dnase_ideas_df.shape[0]:
    print "Great! dnase ideas df are perfectly concatenated"
else:
    warnings.warn("Please, re-check the concatenation...")

""" Preparing for significant peak and background peaks overlap """
# dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)
len(tf_file_list)

concat_tf_df_list = []
for each_file in tf_file_list:
    #tf_ideas_dict = {}
    tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
    print "\nCurrently Processing %s\n" %(tf_name)
    df_read = pd.read_csv(each_file, sep="\t", header=None)
    df_read = df_read.iloc[:,[0,1,2]]
    df_read.columns = ["tf_chrom", "start", "end"]
    df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
    df_read["tf_start"] = df_read["mid_point"].astype(int) - (50)
    df_read["tf_end"] = df_read["mid_point"].astype(int) + (50)
    df_read["tf_name"] = tf_name

    """ Preparing the file for significant peak overlap"""
    select_cols = ["tf_chrom","tf_start","tf_end","tf_name"]
    df_read = df_read.loc[:,select_cols]
    print df_read

    # append the df of each tf to make the single combined list of bedfiles:
    concat_tf_df_list.append(df_read)

combined_tf_df = pd.concat(concat_tf_df_list, ignore_index=True)
combined_tf_sorted_df = combined_tf_df.sort_values(["tf_chrom", "tf_start"])


""" Intersection b/w dnase_ideas_pybed and combined_tf_list """
tf_bedfile = pybedtools.BedTool.from_dataframe(combined_tf_sorted_df)
print "\nPreparing the intersection b/w combined tf df and ideas each sites and peaks"

pybed_outfile = join(output_dir, ("pybed_tf_merged_ideas_intersect_test.bed"))
tf_ideas_intersect = tf_bedfile.intersect(dnase_ideas_pybed, wa=True, wb=True, output=pybed_outfile)
print tf_ideas_intersect.head()

tf_ideas_intersect_df = pd.read_csv(tf_ideas_intersect.fn, sep="\t", header=None)
tf_ideas_intersect_df = tf_ideas_intersect_df.iloc[:,0:9]
tf_ideas_intersect_df.head()
tf_ideas_intersect_df.columns = ["tf_chrom","tf_start","tf_end","tf_name", "ideas_chrom", "ideas_start", "ideas_end", "ideas_state", "state_anno"]
tf_ideas_grouped_df = tf_ideas_intersect_df.groupby(["ideas_chrom", "ideas_start", "ideas_end", "ideas_state", "state_anno", "tf_name"]).count()

tf_ideas_grouped_df.to_csv(join(output_dir,"final_tf_hotspots_sites_with_TF_details_others_for_dnase.bed"), sep="\t", header=True, index=True)
with open(join(output_dir,"final_tf_hotspots_sites_with_TF_details.bed.pkl"), "w") as outfile:
    pickle.dump(tf_ideas_grouped_df, outfile)

df_selected = tf_ideas_grouped_df.reset_index().iloc[:,0:6]
#df_selected["tf_counts"] = df_selected.groupby(["ideas_chrom", "ideas_start", "ideas_end"])["tf_name"].transform("count") ## gives counts of TF as well
df_selected.to_csv(join(output_dir,"final_tf_hotspots_sites_with_TF_details_for_dnase.bed"), sep="\t", header=True, index=True)
df_selected_grouped = df_selected.groupby(["ideas_chrom", "ideas_start", "ideas_end", "ideas_state","state_anno"]).count()
df_ordered_by_tf_count = df_selected_grouped.reset_index().rename(columns={"tf_name" : "tf_counts"})
df_sorted_by_tf_count = df_ordered_by_tf_count.sort_values("tf_counts", ascending=False)
tf_bound_site_distribution = df_sorted_by_tf_count["ideas_state"].value_counts()

tf_bound_site_distribution.to_csv(join(output_dir,"final_dnase_tf_associated_piechart_distribution_data.txt"), sep="\t", header=True, index=True)
df_sorted_by_tf_count.to_csv(join(output_dir,"final_tf_hotspots_sites_distribution_sorted_for_dnase.bed"), sep="\t", header=True, index=True)
with open(join(output_dir,"final_tf_hotspots_sites_distribution.pkl"), "w") as outfile:
    pickle.dump(df_ordered_by_tf_count, outfile)

""" Details for 1 or more bound sites """
#df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
df_sorted_by_tf_count["segment_size"] = (df_sorted_by_tf_count["ideas_end"] - df_sorted_by_tf_count["ideas_start"]) + 1
df_sorted_by_tf_count["counter"] = 1
df_sorted_by_tf_count_gp = df_sorted_by_tf_count.groupby(["tf_counts"]).apply(lambda x: (x["segment_size"].sum(), x["counter"].sum()))
df_tf_count_reset = df_sorted_by_tf_count_gp.reset_index()
df_tf_count_reset[["segment_size_sum", "total_sites"]] = df_tf_count_reset.loc[:,0].apply(pd.Series)
df_tf_count_reset["tf_counts_per_kb_segment"] = (df_tf_count_reset["tf_counts"]*df_tf_count_reset["total_sites"]/df_tf_count_reset["segment_size_sum"]) * 1000 #tf_per_kb = (tf_counts/seg_size) * 1000
select_cols = ["tf_counts", "total_sites", "segment_size_sum", "tf_counts_per_kb_segment"]
final_tf_count_intersect_df = df_tf_count_reset.loc[:, select_cols]

""" Find the details for zero TF sites """
pybed_outfile_v = join(output_dir, ("pybed_tf_merged_ideas_outersect_test.bed"))
tf_ideas_outersect = dnase_ideas_pybed.intersect(tf_bedfile, wa=True, wb=True, v=True, output=pybed_outfile_v)
print tf_ideas_outersect.head()
df_outersect = pd.read_csv(tf_ideas_outersect.fn, sep="\t", header=None)
df_outersect = df_outersect.iloc[:,0:5]
df_outersect.columns = ["chrom", "start", "end", "state", "anno_state"]
df_outersect.head()
df_outersect.shape[0]

# needed for R to fit in the data frame for the barplot:
prom_assoc_count = df_outersect[df_outersect["anno_state"] == "prom_assoc"].shape[0]; prom_assoc_count
sEnh_assoc_count = df_outersect[df_outersect["anno_state"] == "sEnh_assoc"].shape[0]; sEnh_assoc_count
wEnh_assoc_count = df_outersect[df_outersect["anno_state"] == "wEnh_assoc"].shape[0]; wEnh_assoc_count

df_outersect["segment_size"] = (df_outersect["end"] - df_outersect["start"]) + 1
df_outersect_total_size = sum(df_outersect["segment_size"])
df_outersect_norm = pd.DataFrame({"tf_counts": [0], "total_sites": [df_outersect.shape[0]], "segment_size_sum" : [df_outersect_total_size]})
df_outersect_norm["tf_counts_per_kb_segment"] = (df_outersect_norm["tf_counts"]*df_outersect_norm["total_sites"]/df_outersect_norm["segment_size_sum"]) * 1000 #tf_per_kb = (tf_counts/seg_size) * 1000
select_cols = ["tf_counts", "total_sites", "segment_size_sum", "tf_counts_per_kb_segment"]
final_norm_outersect_df = df_outersect_norm.loc[:, select_cols]

""" Merge the zero Tf sites and 1 or more bound sites for final data """
combined_tf_count_distribution =  pd.concat([final_norm_outersect_df, final_tf_count_intersect_df], ignore_index=True)
combined_tf_count_distribution.to_csv(join(output_dir,"final_tf_hotspots_double_panels_barplot_data.bed"), sep="\t", header=True, index=True)

""" Merge the zero TF and 1 or more bound sites for typical single panel barplot_data """
dnase_tf_bound = df_sorted_by_tf_count.groupby(["tf_counts", "state_anno"]).size().reset_index()
dnase_tf_bound.rename(columns={0:"total_sites"}, inplace=True)
dnase_tf_unbound = pd.DataFrame({"tf_counts": [0,0,0] , 
                            "state_anno":["prom_assoc","sEnh_assoc","wEnh_assoc"], 
                            "total_sites": [prom_assoc_count, sEnh_assoc_count, wEnh_assoc_count]
                            })
dnase_bound_unbound = pd.concat([dnase_tf_unbound, dnase_tf_bound], ignore_index=True)
dnase_bound_unbound.to_csv(join(output_dir,"final_dnase_tf_hotspots_single_barplot_data.bed"), sep="\t", header=True, index=True)



