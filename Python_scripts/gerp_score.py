import pandas as pd
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from pybedtools import BedTool
from os.path import join
from os.path import basename
from os.path import splitext
import warnings


ideas_tf_hotspots_file = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/files/genomic_sites/final_tf_hotspots_sites_distribution_100bp_merged_with_coordinates.txt")
gerp_chrom_wise_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/files/hg19_GERP_scores_perbase/*.gerp")
gerp_chrom_wise_files = glob.glob(gerp_chrom_wise_dir)

output_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/output")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


""" Select the limits for high TF bound regions"""
lower_limit=50
upper_limit=125

hotspot_file = pd.read_csv(ideas_tf_hotspots_file, sep="\t", index_col=0)
hotspot_df = hotspot_file.loc[(hotspot_file["tf_counts"] >=lower_limit) & (hotspot_file["tf_counts"] <=upper_limit),:]

gerp_df_list = []

gerp_whole_genome_file = join(output_dir, "hg19_GERP_scores_per_base_whole_genome.txt")
if not.os.path.exists(gerp_whole_genome_file):
    for each_chrom_file in gerp_chrom_wise_files:
        gerp_df = pd.read_csv(each_chrom_file, sep="\t", header=None)
        #gerp_chunk_df = pd.read_csv(gerp_chrom_wise_files[1], sep="\t", header=None, chunksize=500000)
        #gerp_df = pd.concat(gerp_chunk_df, ignore_index=True)
        gerp_df.columns = ["chrom_start", "neutral_rate", "rs_scores"]

        chrom_name = basename(each_chrom_file).split(".")[0]
        print "processing chromosome ...", chrom_name
        gerp_df["chrom"] = chrom_name
        gerp_df["chrom_end"] = gerp_df["chrom_start"] + 1
        gerp_df = gerp_df.loc[:, ["chrom", "chrom_start","chrom_end", "neutral_rate", "rs_scores"]]
        gerp_df.to_csv(join(output_dir, "hg19_GERP_scores_per_base_" + chrom_name + ".txt"), sep="\t", header=True, index=False)
        gerp_df_list.append(gerp_df)

    gerp_whole_genome_df = pd.concat(gerp_df_list, ignore_index=True)
    #gerp_whole_genome_sorted_df = gerp_whole_genome_df.sort_values(["chrom"]).reset_index(drop=True)
    gerp_whole_genome_df.to_csv(join(output_dir, "hg19_GERP_scores_per_base_whole_genome.txt"), sep="\t", header=True, index=False)
    print "hg19_GERP_scores_per_base file generation successfully completed !!!"

else:
    print "hg19_GERP_scores_per_base file already exists"

gerp_file = BedTool("hg19_GERP_scores_per_base_whole_genome.txt")
hotspot_file = BedTool("TF_hotspot_sites.txt")
hotspot_gerp = hotspot_file.intersect(gerp_file, wa =True, wb = True, output = "TF_hotspot_hg19_gerp_score_intersect.txt")
hotspot_gerp_groupby = hotspot_gerp.groupby([1,2,3], c=11, o=["mean"])
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


