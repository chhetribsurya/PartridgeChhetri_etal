import pandas as pd
import pybedtools
from pybedtools import BedTool
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

#output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/assign_IDEAS_states"
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/assign_IDEAS_states/ryne_data"
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
dnase_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/assign_IDEAS_states/ryne_data"
#dnase_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/assign_IDEAS_states"
#dnase_rep1_bed_file = join(dnase_file_dir, "chris_data.txt")
dnase_rep1_bed_file = join(dnase_file_dir, "Great_analysis_ryne_hotspot.txt")

dnase_rep1_df = pd.read_csv(dnase_rep1_bed_file, sep = "\t", header = None)
dnase_rep1_df = dnase_rep1_df.iloc[:, [0, 1, 2]]
dnase_rep1_df.columns = ["chrom", "start", "end"]

### Intersect Hepg2 bed file and Filter the rows with max base overlap size for duplicated rows:
hepg2_ideas_df = pd.read_csv(ideas_hepg2_file, sep="\t")
hepg2_ideas_df = hepg2_ideas_df.iloc[:, 1:5]
hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "chromStart"])

hepg2_bed_file = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
dnase_bed_file = pybedtools.BedTool.from_dataframe(dnase_rep1_df)

dnase_ideas_intersect = dnase_bed_file.intersect(hepg2_bed_file, wao=True)

### remove the feature with 0 overlap size(last column)
dnase_ideas_intersect_df = pd.read_csv(dnase_ideas_intersect.fn, sep="\t", header=None)
intersect_df =  dnase_ideas_intersect_df[dnase_ideas_intersect_df.iloc[:, -1] > 0] 

### Filter the rows with max base overlap size for duplicated rows:(column 8 here represents the size for base overlap)
duplicated_df = intersect_df[intersect_df.duplicated([0,1,2], keep=False)]
duplicated_filtered_df = duplicated_df.loc[duplicated_df.groupby([0,1,2])[7].idxmax()]
non_duplicated_df = intersect_df[~intersect_df.duplicated([0,1,2], keep=False)] # non-duplicated df
if (intersect_df.shape[0] == duplicated_df.shape[0] + non_duplicated_df.shape[0]): # 47743 = 34578 + 13165 
    print "Duplicated Rows filtered after groupby and merged successfully"

dnase_uniq_df = pd.concat([duplicated_filtered_df, non_duplicated_df], ignore_index=True)
dnase_final_df = dnase_uniq_df.iloc[:,range(0,3)+[6]]
dnase_final_df.columns = ["chrom", "start", "end", "state"]
dnase_final_df = dnase_final_df.sort_values(["chrom", "start"])
dnase_final_df = dnase_final_df.reset_index(drop=True) # reordering of the index
dnase_ideas_distribution = dnase_final_df["state"].value_counts()

dnase_ideas_distribution.to_csv(join(output_dir, "uniq_peaks_wholegenome_piechart_distribution_data.txt"), header=True, index=True, sep="\t")
dnase_final_df.to_csv(join(output_dir, "uniq_peaks_annotation_with_ideas_states.bed"), header=True, index=True, sep="\t")


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
    dnase_ideas_sorted_df.to_csv(join(output_dir, Merged_state_anno + "_peaks_ideas_states.bed"), header=True, sep="\t")

    """ Generate pybed object"""
    dnase_ideas_pybed = pybedtools.BedTool.from_dataframe(dnase_ideas_sorted_df) 
    #merge_ideas_pybed = dnase_ideas_pybed.merge(d=1, c=[4,5], o=["distinct", "distinct"])
    return(dnase_ideas_pybed)

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
    pybed_df.columns = ["chrom", "start","end", "state", "state_anno", "strand"]
    concat_list_df.append(pybed_df)
    df_line_count.append(pybed_df.shape[0])

dnase_ideas_df = pd.concat(concat_list_df, ignore_index=True)
dnase_ideas_df = dnase_ideas_df.sort_values(["chrom", "start"]).reset_index(drop=True)
dnase_ideas_df_final = dnase_ideas_df.iloc[:, range(0,5)]
dnase_ideas_df_final.to_csv(join(output_dir, "uniq_peaks_with_prom_sEnh_wEnh_category.txt"), sep="\t", header=True, index=True)
