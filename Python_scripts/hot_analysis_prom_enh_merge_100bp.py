import pandas as pd
import pybedtools
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from rpy2.robjects.packages import importr  #stats = importr('stats')
from rpy2.robjects.vectors import FloatVector
from os.path import join
from os.path import basename
from os.path import splitext


###################################
###################################

""" Read hepg2 ideas whole genome file """
hepg2_df = pd.read_csv("/gpfs/gpfs1/home/schhetri/for_chris/HepG2_ideas_whole_genome.bed", sep="\t")
hepg2_df.columns = ["chrom", "start", "end", "state", "signal", "strand" ] 
hepg2_df

sub_output_dir = "/gpfs/gpfs1/home/schhetri/for_chris"
if not os.path.exists(sub_output_dir):
	os.makedirs(sub_output_dir)

""" HepG2 Ideas Whole genome sorted"""
hepg2_df_sorted = hepg2_df.sort_values(["chrom", "start"])
hepg2_df_sorted.to_csv(join(sub_output_dir, "HepG2_ideas_whole_genome_sorted.bed"), sep="\t", header=False, index=False)

""" HepG2 Ideas Promoter and enhancer associated region sorted """
prom_enh_df_sorted = hepg2_df[hepg2_df["state"].isin(target_state)].sort_values(["chrom","start"])
prom_enh_df_sorted.to_csv(join(sub_output_dir, "HepG2_ideas_prom_enh_sorted.bed"), sep="\t", header=False, index=False)


###################################
###################################


# ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

""" Select the state number, to be analysed on """
select_state_num = range(1,8+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
#select_state_num = range(1,9+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
#select_state_num = range(1,36+1) # Add +1; range() is exclusive:

""" Generate file list"""
ideas_file_dir = "/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state"
hepg2_file_list = [ each_file for each_file in glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) 
				if os.path.basename(each_file).startswith("HepG2")]

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_prom_enh/merged_d_100"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

""" Read the mnemonics file, and create a dictionary of states with state number ID"""
read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
state_list = read_file["Mnemonics"]
state_num_list = [ i for i in range(1,len(state_list)+1) ]
ideas_state_dict = dict(zip(state_num_list,state_list))
target_state = [ideas_state_dict[each] for each in select_state_num] #Depends on the select_state_num variable input


def merge_ideas_pybed(Ideas_mnemonics_file, File_to_merge_list, Select_state_num):
	ideas_mnemonics_file = Ideas_mnemonics_file
	file_list = File_to_merge_list
	select_state_num = Select_state_num

	""" Depends on the select_state_num variable input"""
	target_state = [ideas_state_dict[each] for each in select_state_num]
	# print target_state
	ideas_select_df_list = []
	for each_ideas_file in file_list:
		cell_line_name = basename(each_ideas_file).split("_")[0]	
		ideas_df = pd.read_csv(each_ideas_file, sep="\t")
		ideas_df = ideas_df.iloc[:,[1,2,3,4]]

		ideas_df.columns = ["chrom", "start", "end", "state"]	
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
hepg2_ideas_pybed =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, select_state_num) 
# df = pd.read_csv(hepg2_file_list[0], sep="\t") ; # In [779]: df.shape; Out[779]: (2622465, 10)


""" Preparing for significant peak and background peaks overlap """
# dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)


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


""" Intersection b/w hepg2_ideas_pybed and combined_tf_list """
tf_bedfile = pybedtools.BedTool.from_dataframe(combined_tf_sorted_df)
print "\nPreparing the intersection b/w combined tf df and ideas each sites and peaks"

pybed_outfile = join(output_dir, ("final_tf_merged_ideas_intersect_100bp.bed"))
tf_ideas_intersect = tf_bedfile.intersect(hepg2_ideas_pybed, wa=True, wb=True, output=pybed_outfile)
print tf_ideas_intersect.head()

tf_ideas_intersect_df = pd.read_csv(tf_ideas_intersect.fn, sep="\t", header=None)
tf_ideas_intersect_df = tf_ideas_intersect_df.iloc[:,0:8]
tf_ideas_intersect_df.columns = ["tf_chrom","tf_start","tf_end","tf_name", "ideas_chrom", "ideas_start", "ideas_end", "ideas_state"]
tf_ideas_grouped_df = tf_ideas_intersect_df.groupby(["ideas_chrom", "ideas_start", "ideas_end", "ideas_state", "tf_name"]).count()

tf_ideas_grouped_df.to_csv(join(output_dir,"final_tf_hotspots_sites_with_TF_details_others_100_bp.bed"), sep="\t", header=True, index=True)
with open(join(output_dir,"final_tf_hotspots_sites_with_TF_details.bed.pkl"), "w") as outfile:
	pickle.dump(tf_ideas_grouped_df, outfile)

df_selected = tf_ideas_grouped_df.reset_index().iloc[:,0:5]
df_selected.to_csv(join(output_dir,"final_tf_hotspots_sites_with_TF_details_100bp.bed"), sep="\t", header=True, index=True)
df_selected_grouped = df_selected.groupby(["ideas_chrom", "ideas_start", "ideas_end", "ideas_state"]).count()
df_ordered_by_tf_count = df_selected_grouped.reset_index().rename(columns={"tf_name" : "tf_counts"})
df_sorted_by_tf_count = df_ordered_by_tf_count.sort_values("tf_counts", ascending=False)

df_ordered_by_tf_count.to_csv(join(output_dir,"final_tf_hotspots_sites_distribution_100bp.bed"), sep="\t", header=True, index=True)
df_sorted_by_tf_count.to_csv(join(output_dir,"final_tf_hotspots_sites_distribution_sorted_100bp.bed"), sep="\t", header=True, index=True)
with open(join(output_dir,"final_tf_hotspots_sites_distribution.pkl"), "w") as outfile:
	pickle.dump(df_ordered_by_tf_count, outfile)

""" Details for 1 or more bound sites """
#df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
df_sorted_by_tf_count["segment_size"] = (df_sorted_by_tf_count["ideas_end"] - df_sorted_by_tf_count["ideas_start"]) + 1
df_sorted_by_tf_count["counter"] = 1
df_sorted_by_tf_count_gp = df_sorted_by_tf_count.groupby(["tf_counts"]).apply(lambda x: (x["segment_size"].sum(), x["counter"].sum()))
df_tf_count_reset = df_sorted_by_tf_count_gp.reset_index()
df_tf_count_reset[["segment_size_sum", "total_sites"]] = df_tf_count_reset.loc[:,0].apply(pd.Series)
df_tf_count_reset["tf_counts_per_kb_segment"] = df_tf_count_reset["tf_counts"]/df_tf_count_reset["segment_size_sum"] * 1000 #tf_per_kb = (tf_counts/seg_size) * 1000
select_cols = ["tf_counts", "total_sites", "segment_size_sum", "tf_counts_per_kb_segment"]
final_tf_count_intersect_df = df_tf_count_reset.loc[:, select_cols]

""" Find the details for zero TF sites """
pybed_outfile_v = join(output_dir, ("final_tf_merged_ideas_outersect_100.bed"))
tf_ideas_outersect = hepg2_ideas_pybed.intersect(tf_bedfile, wa=True, wb=True, v=True, output=pybed_outfile_v)
print tf_ideas_outersect.head()
df_outersect = pd.read_csv(tf_ideas_outersect.fn, sep="\t", header=None)
df_outersect = df_outersect.iloc[:,0:5]
df_outersect.columns = ["chrom", "start", "end", "state", "cellline_state"]

df_outersect["segment_size"] = (df_outersect["end"] - df_outersect["start"]) + 1
df_outersect_total_size = sum(df_outersect["segment_size"])
df_outersect_norm = pd.DataFrame({"tf_counts": [0], "total_sites": [df_outersect.shape[0]], "segment_size_sum" : [df_outersect_total_size]})
df_outersect_norm["tf_counts_per_kb_segment"] = df_outersect_norm["tf_counts"]/df_outersect_norm["segment_size_sum"] * 1000 #tf_per_kb = (tf_counts/seg_size) * 1000
select_cols = ["tf_counts", "total_sites", "segment_size_sum", "tf_counts_per_kb_segment"]
final_norm_outersect_df = df_outersect_norm.loc[:, select_cols]

""" Merge the zero Tf sites and 1 or more bound sites for final data """
combined_tf_count_distribution =  pd.concat([final_norm_outersect_df, final_tf_count_intersect_df], ignore_index=True)
combined_tf_count_distribution.to_csv(join(output_dir,"final_tf_hotspots_barplot_data_100bp.bed"), sep="\t", header=True, index=True)

#####################

""" For the details on zero TF binding site distribution""" 

""" Total sites for each state of prom and enh associated """
df_ori = pd.read_csv(hepg2_ideas_pybed.fn, sep="\t", header=None)
df_ori.columns = ["chrom", "start", "end", "state", "cellline_state"]

total_count_dict = {}
total_coord_dict = {}
for each in target_state:
	total_count_df = df_ori[df_ori["state"].apply( lambda x : each in x.split(","))]
	total_coord_dict[each] = total_count_df
	total_count_dict[each] = total_count_df.shape[0]


""" For unbound sites """ # Find the non intersection:
df_outersect = pd.read_csv(tf_ideas_outersect.fn, sep="\t", header=None)
df_outersect.columns = ["chrom", "start", "end", "state", "cellline_state"]
zero_tf_count = df_outersect.shape[0]
zero_tf_count

unbound_dict = {}
unbound_dict_percent = {}
unbound_coord_dict = {}
for each in target_state:
	outersect_df = df_outersect[df_outersect["state"].apply( lambda x : each in x.split(","))]
	unbound_coord_dict[each] =  outersect_df
	unbound_dict[each] = outersect_df.shape[0]
	unbound_dict_percent[each] = (outersect_df.shape[0]/float(total_count_dict[each])) * 100


""" For bound sites """ # use the final df produced from the analysis:
df_sorted_by_tf_count.columns = ["chrom", "start", "end", "state", "tf_counts", "seg_size", "counter"]

bound_dict = {}
bound_dict_percent = {}
bound_coord_dict = {}
for each in target_state:
	intersect_df = df_sorted_by_tf_count[df_sorted_by_tf_count["state"].apply( lambda x : each in x.split(","))]
	bound_coord_dict[each] =  intersect_df
	bound_dict[each] = intersect_df.shape[0]
	bound_dict_percent[each] = (intersect_df.shape[0]/float(total_count_dict[each])) * 100


""" Forming the dataframe using the above values """
final_bound_unbound_df = pd.DataFrame({"state": unbound_dict.keys(), "total_sites": total_count_dict.values(), "bound_sites": bound_dict.values(), 
			"unbound_sites": unbound_dict.values(), "bound_percent": bound_dict_percent.values(), "unbound_percent": unbound_dict_percent.values() })
final_bound_unbound_df = final_bound_unbound_df.loc[:, ["state", "total_sites", "bound_sites", "unbound_sites", "bound_percent", "unbound_percent"]]
final_bound_unbound_df
final_bound_unbound_df.to_csv(join(output_dir, "final_bound_unbound_site_distribution_100bp.bed"), header=True, sep="\t")

# just a double confirmation of those total sites and total percent contribution:
final_bound_unbound_df["total"] = final_bound_unbound_df["bound_sites"] + final_bound_unbound_df["unbound_sites"]
final_bound_unbound_df["total_percent"] = final_bound_unbound_df["bound_percent"] + final_bound_unbound_df["unbound_percent"]
final_bound_unbound_df

for each in target_state:
    total_sites = bound_dict[each] + unbound_dict[each]
    print each, ":", total_sites

""" But for merged ones, need to calculate the following way """
total_hepg2_sites = hepg2_ideas_pybed.count(); total_hepg2_sites
total_tf_covered_sites = df_sorted_by_tf_count.shape[0]; total_tf_covered_sites
total_0_tf_count_sites = total_hepg2_sites - total_tf_covered_sites; total_0_tf_count_sites
