import pandas as pd
import numpy as np
import seaborn as sns; sns.set()
# %matplotlib # morgan/servers usually get upset by the backends
import os
import re
import collections
import pybedtools
import glob
import json
import pickle
from os.path import splitext
from os.path import basename	
from os.path import join

###################################

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_co_binding_total"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

""" Choose the file category or TF category you want to analyse on """
#file_category = "all_tf" #file_category = "all_tf"; file_category = "cr_tf"; file_category = "dbf_tf"

""" All TF containing dir """
#dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob.glob(dir_path) # represents all tf file list

###################################

"""Read multiple excel sheets from the same file"""
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx"
df_xls = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")

df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
xls_tf_list =  df_xls["Target"].tolist()

cr_df = df_xls[df_xls["Category"] == "CR/CF"]
dbf_df = df_xls[~(df_xls["Category"] == "CR/CF")]

""" Check the xls TF list and the file suffix tf_list to make sure, if they are in the TF list"""
TF_check_list = []
for each_file in all_tf_file_list:
	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
	TF_check_list.append(tf_name)

# df_xls[df_xls["Target"].isin(TF_check_list)] # or,
if sorted(TF_check_list) == sorted(xls_tf_list):
	print "\nGreat!! TF list resembles b/w the file list and xls tf list..."
else:
	print "\nWarning: Your files TF list doesn't resemble the xls tf list..."


""" Select the category of TF for chosing the file list """
cr_tf_file_list = []
for each_tf in cr_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			cr_tf_file_list.append(each_file)

dbf_tf_file_list = []
for each_tf in dbf_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			dbf_tf_file_list.append(each_file)

def combine_tf_bedfiles_midpt_centered(TF_file_list_of_interest):
	concat_file_list = []
	for each_file in tf_file_list:
		tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
		df_read = pd.read_csv(each_file, header=None, sep="\t")
		df_read = df_read.iloc[:,[0,1,2]]
		df_read.columns = ["peak_chrom", "start", "end"]
		df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
		df_read["peak_start"] = df_read["mid_point"].astype(int)
		df_read["peak_end"] = df_read["mid_point"].astype(int) + 1
		df_read["strand"] = "."
		df_read["tf_name"] = tf_name
		select_cols = ["peak_chrom","peak_start","peak_end","start","end", "strand", "tf_name"]
		df_read = df_read.loc[:,select_cols]
		print "test\n\n"
		print df_read
		concat_file_list.append(df_read)

	combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
	combined_tf_df.to_csv(join(output_dir, "final_combined_tf_coord_df.bed"), header=False, index=False, sep="\t")
	return(combined_tf_df)

# combined_tf_df = combine_tf_bedfiles_midpt_centered(tf_file_list)


def combine_tf_bedfiles(TF_file_list_of_interest):
	tf_file_list = TF_file_list_of_interest
	concat_file_list = []
	for each_file in tf_file_list:
		tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
		df_read = pd.read_csv(each_file, header=None, sep="\t")
		df_read = df_read.iloc[:,[0,1,2]]
		df_read.columns = ["peak_chrom", "start", "end"]
		df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
		df_read["peak_start"] = df_read["mid_point"].astype(int) - (50)
		df_read["peak_end"] = df_read["mid_point"].astype(int) + (50)
		df_read["strand"] = "."
		df_read["tf_name"] = tf_name
		select_cols = ["peak_chrom","peak_start","peak_end","start","end", "strand", "tf_name"]
		df_read = df_read.loc[:,select_cols]
		print "test\n\n"
		print df_read
		concat_file_list.append(df_read)

	combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
	combined_tf_df.to_csv(join(output_dir, "final_combined_tf_coord_df.bed"), header=False, index=False, sep="\t")
	return(combined_tf_df)

#combined_tf_df = combine_tf_bedfiles(all_tf_file_list)


def load_combined_peaks_coord_pybedtool_object(file_name_with_full_path): 
	each_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["each_file"] = each_file
	sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
	os.environ["sorted_peak_file"] = sorted_peak_file
	# if not os.path.exists(join(output_dir, sorted_peak_file)):
	CMD = 'sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file'
	os.system(CMD)   
  
	print "Generating bedtools object..."
	peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
	return(peak_bed)

#combined_tf_bedfile = load_combined_peaks_coord_pybedtool_object(join(output_dir, "final_combined_tf_coord_df.bed"))
# pybed_df = pd.read_csv(combined_tf_bedfile.fn, sep="\t")



def cobinding_tf_peaks_count_perc_overlap_info(TF_file_list_of_interest, combined_TF_bedfile):
	tf_file_list = TF_file_list_of_interest
	combined_tf_bedfile = combined_TF_bedfile
	master_tf_dict = {}
	master_check_list = []
	master_name = []

	for each_file in tf_file_list:
		tf_cobind_dict = {}
		tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
		print "\nCurrently Processing %s\n" %(tf_name)
		df_read = pd.read_csv(each_file, sep="\t", header=None)
		df_read = df_read.iloc[:,[0,1,2]]
		df_read.columns = ["tf_chrom", "start", "end"]

		df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
		df_read["tf_start"] = df_read["mid_point"].astype(int) - (50)
		df_read["tf_end"] = df_read["mid_point"].astype(int) + (50)
		df_read["tf_name"] = tf_name
		select_cols = ["tf_chrom","tf_start","tf_end","tf_name"]
		df_read = df_read.loc[:,select_cols]

		tf_bedfile = pybedtools.BedTool.from_dataframe(df_read)
		print tf_bedfile.head()

		"""Pybedtool intersection of tf_bedfile and combined_tf_bedfile"""
		print "Processing the intersection for", tf_name
		pybed_outfile = join(output_dir, (tf_name + "_and_combined_tf_intersect.bed"))
		tf_tf_intersect = tf_bedfile.intersect(combined_tf_bedfile, wa=True, wb=True, output=pybed_outfile)
		print tf_tf_intersect.head()

		pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
		header = ["tf_chrom", "tf_start", "tf_end", "tf_name", "peak_chrom", "peak_start", "peak_end", "start", "end", "strand", "cobind_tf_name"]
		pybed_df.columns = header

		"""Comment pybed_df if you need the original TF or the self_Tf-self_Tf intersect too"""
		pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] 
		cobinding_tf_list = pybed_df.loc[:,"cobind_tf_name"].unique().tolist()

		# Add itself to the list of complex, but comment the line if you need self_TF-self_TF intersect
		cobinding_tf_list.append(tf_name) # since; pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] commented;

		# Find unique non-repeating peaks wrt to chrom, start, tf_name, cobind_tf_name:  
		pybed_df_group = pybed_df.groupby(["tf_chrom","tf_start","tf_end","tf_name","cobind_tf_name"]).count().reset_index()

		# Groupby tf_name, and cobind_tf_name to get total peak_overlap_count:
		cobinding_tf_peaks_count = pybed_df_group.groupby(["tf_name","cobind_tf_name"]).count().iloc[:,0].reset_index()
		cobinding_tf_peaks_count = cobinding_tf_peaks_count.rename(columns={"tf_chrom" : "peak_overlap_count"})
		cobinding_tf_peaks_count["total_peak_count"] = df_read.shape[0] # initial or original peak count
		cobinding_tf_peaks_count["percent_overlap"] = (cobinding_tf_peaks_count["peak_overlap_count"]/cobinding_tf_peaks_count["total_peak_count"])*100
		cobinding_tf_peaks_count = cobinding_tf_peaks_count.sort_values(["percent_overlap"], ascending=False)

		tf_cobind_dict["cobind_complex"] = pybed_df
		tf_cobind_dict["cobind_list"] = cobinding_tf_list
		tf_cobind_dict["cobind_peak_count"] = cobinding_tf_peaks_count

		master_tf_dict[tf_name] = tf_cobind_dict
		master_check_list.append(cobinding_tf_list)
		master_name.append(tf_name)

	return(master_tf_dict)


""" Final Cobind peak count dataframe """
def cobinding_tf_peaks_overlap_final_output(master_tf_dict_info, File_category):
	master_TF_dict = master_tf_dict_info
	file_suffix =  File_category
	cobind_peakcount_df_list = []
	for each_tf in master_TF_dict:
		cobind_peakcount_df_list.append(master_TF_dict[each_tf]["cobind_peak_count"])

	cobind_peakcount_dframe = pd.concat(cobind_peakcount_df_list)
	cobind_peak_count_df = cobind_peakcount_dframe.set_index(["tf_name", "cobind_tf_name"]) ### For heirarchial multiple indexing
	cobind_peak_count_df.columns = ["overlap_count", "total_peaks", "percent_overlap"]
	cobind_peak_count_df.to_csv(join(output_dir,"final_cobinding_peak_count_perc_overlap_" + file_suffix + ".bed"), sep="\t", header=True, index=True)

	cobind_peak_count_sorted = cobind_peak_count_df.sort_values(["percent_overlap"], ascending = False)
	cobind_peak_count_sorted.to_csv(join(output_dir,"final_cobinding_peak_count_perc_overlap_sorted_" + file_suffix + ".bed"), sep="\t", header=True, index=True)

	### Generating a heatmap data from the cobind peakcount dataframe:
	peak_table = cobind_peakcount_dframe.iloc[:,[0,1,4]]
	peak_heatmap_table = peak_table.pivot(index="tf_name", columns="cobind_tf_name", values = "percent_overlap")
	peak_heatmap_table.to_csv(join(output_dir,"final_cobinding_peak_count_heatmap_" + file_suffix + ".bed"), sep="\t", header=True, index=True)

	return(cobind_peak_count_df,cobind_peak_count_sorted)


# def main()
""" Choose the file category or TF category you want to analyse on """
file_category = "dbf_tf" # cr_tf; dbf_tf, all_tf
#file_category = "cr_tf" # cr_tf; dbf_tf, all_tf
#file_category = "all_tf" # cr_tf; dbf_tf, all_tf
output_dir = join(main_dir, file_category)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

"""Note: Make sure that tf_file_list matches with the file_category""" 
"""For "test_tf" & "all_tf" file category, all_tf_file_list works or is common to both"""
combined_tf_df = combine_tf_bedfiles(dbf_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list
#combined_tf_df = combine_tf_bedfiles(cr_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list
#combined_tf_df = combine_tf_bedfiles(all_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list
combined_tf_bedfile = load_combined_peaks_coord_pybedtool_object(join(output_dir, "final_combined_tf_coord_df.bed"))

master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(dbf_tf_file_list, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list
#master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(cr_tf_file_list, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list
#master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(all_tf_file_list, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list

final_cobind_df,final_cobind_df_sorted = cobinding_tf_peaks_overlap_final_output(master_TF_dict, file_category) 

print "Job completed successfully!!"





