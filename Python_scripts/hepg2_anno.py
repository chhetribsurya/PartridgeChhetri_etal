import pandas as pd
import pybedtools
import pickle
import scipy 
import os
import re
import glob
from os.path import join
from os.path import basename
from os.path import splitext
from os.path import expanduser
from 
import warnings

#ideas_file_dir = "/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state"
ideas_file_dir = expanduser("~/Dropbox/for_chris")
ideas_file = join(ideas_file_dir, "hepg2_ideas_36_dense.bed")

#output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_prom_enh_categorised/merged_d_100"
output_dir = join(ideas_file_dir, "ideas_annotation")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

read_ideas_df = pd.read_csv(ideas_file, sep="\t", header=None)
ideas_df_select = read_ideas_df.iloc[:,[0,1,2,3]]
ideas_df_select.columns = ["chrom", "start", "end", "ideas_state"]


### For whole genome annotation:
map_ideas = {"Tss" : "Prom assoc", "TssF" :  "Prom assoc", 'TssCtcf': "Prom assoc", 'PromCtcf': "Prom assoc", 
"TssW" : "Prom assoc", "PromP": "Prom assoc", 'PromF1': "Prom assoc", 'PromF2': "Prom assoc",
"Pol2" : "Gene body", "Gen5" : "Gene body", "Gen3": "Gene body", "Gen3Ctcf" : "Gene body", "Elon" : "Gene body", "ElonW" : "Gene body",
"CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", 
"Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh",
"Repr1" : "Heterochrom_repressed", "Repr2" : "Heterochrom_repressed", "ReprD":  "Heterochrom_repressed", "LowReprW":  "Heterochrom_repressed", 
"Low1":  "Heterochrom_repressed", "Low2":  "Heterochrom_repressed", "Quies" :  "Heterochrom_repressed", "Art1":  "Heterochrom_repressed",
"Art2":  "Heterochrom_repressed", "Zero" : "Heterochrom_repressed", 
"DnaseD1": "Euchromatin assoc", "DnaseD2": "Euchromatin assoc", "FaireW1": "Euchromatin assoc", "FaireW2": "Euchromatin assoc"}
ideas_df_select["ideas_anno"] =  ideas_df_select["ideas_state"].map(map_ideas)
ideas_df_select = ideas_df_select.loc[:, ["chrom", "start", "end", "ideas_anno"]]

### For Cis-regulatory elements:
#map_ideas = {"Tss" : "Strong Prom", "TssF" :  "Strong Prom", "TssW" : "Weak Prom", "PromP": "Weak Prom",
#"Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh" } 
#ideas_df_select["ideas_anno"] =  ideas_df_select["ideas_state"].map(map_ideas)

### Selecting only the cpg motif associated to the regulatory regions:
ideas_final = ideas_df_select.dropna(subset=["ideas_anno"])
ideas_final.to_csv(join(output_dir, "Ideas_annotated.bed"), sep = "\t", header = True, index = True )

# ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

""" Select the state number, to be analysed on """
# select_state_num = range(1,8+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
prom_assoc_state_num = range(1,8+1)# Add +1; range() is exclusive:
strong_enh_state_num = [15,16] # Add +1; range() is exclusive:
weak_enh_state_num = range(17,20+1) # Add +1; range() is exclusive:

""" Generate file list"""
ideas_file_dir = "/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state"
hepg2_file_list = [ each_file for each_file in glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) 
				if os.path.basename(each_file).startswith("HepG2")]

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_prom_enh_categorised/merged_d_100"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

""" Read the mnemonics file, and create a dictionary of states with state number ID"""
read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
state_list = read_file["Mnemonics"]
state_num_list = [ i for i in range(1,len(state_list)+1) ]
ideas_state_dict = dict(zip(state_num_list,state_list))

# target_state = [ideas_state_dict[each] for each in select_state_num] #Depends on the select_state_num variable input
prom_assoc_state = [ideas_state_dict[each] for each in prom_assoc_state_num] #Depends on the select_state_num variable input
strong_enh_state = [ideas_state_dict[each] for each in strong_enh_state_num] #Depends on the select_state_num variable input
weak_enh_state = [ideas_state_dict[each] for each in weak_enh_state_num] #Depends on the select_state_num variable input

def merge_ideas_pybed(Ideas_mnemonics_file, File_to_merge_list, Select_state_num, Merged_state_anno):
	ideas_mnemonics_file = Ideas_mnemonics_file
	file_list = File_to_merge_list
	select_state_num = Select_state_num
	merged_state_anno = Merged_state_anno

	""" Depends on the select_state_num variable input"""
	target_state = [ideas_state_dict[each] for each in select_state_num]
	# print target_state
	ideas_select_df_list = []
	for each_ideas_file in hepg2_file_list:
		cell_line_name = basename(each_ideas_file).split("_")[0]	
		ideas_df = pd.read_csv(each_ideas_file, sep="\t")
		ideas_df = ideas_df.iloc[:,[1,2,3,4]]

		ideas_df.columns = ["chrom", "start", "end", "state"]
		ideas_df["state_anno"] = merged_state_anno
		ideas_df["strand"] = "."	
		ideas_select_df = ideas_df[ideas_df["state"].isin(target_state)]	
		ideas_select_df["state_cellLine"] = ideas_select_df["state"] + "_" + cell_line_name
		ideas_select_df_list.append(ideas_select_df)

	combined_ideas_df = pd.concat(ideas_select_df_list, ignore_index=True)
	ideas_sorted_df = combined_ideas_df.sort_values(["chrom", "start"])
	ideas_sorted_df.to_csv(join(output_dir, "final_ideas_states_concat_sorted.bed"), header=True, sep="\t")

	""" Generate pybed object"""
	ideas_pybed = pybedtools.BedTool.from_dataframe(ideas_sorted_df) 
	merge_ideas_pybed = ideas_pybed.merge(d=100, c=[4,5], o=["distinct","distinct"])
	#merge_ideas_pybed = ideas_pybed.merge(d=100, c=[4,5], o=["collapse","collapse"])
	# merge_ideas_df = pd.read_csv(merge_ideas_pybed.fn, sep="\t")
	return(merge_ideas_pybed)

""" Merging the combined_ideas and hepg2 file by calling a function"""
# hepg2_ideas_pybed =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, select_state_num) 
hepg2_ideas_pybed_prom =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, prom_assoc_state_num, "prom_assoc") 
hepg2_ideas_pybed_sEnh =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, strong_enh_state_num, "sEnh_assoc") 
hepg2_ideas_pybed_wEnh =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, weak_enh_state_num, "wEnh_assoc") 

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

hepg2_ideas_df = pd.concat(concat_list_df, ignore_index=True)
hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "start"])
hepg2_ideas_pybed = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
