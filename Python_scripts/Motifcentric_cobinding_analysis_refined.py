import pandas as pd
import numpy as np
import seaborn as sns; sns.set()
import os
import re
import collections
import pybedtools
from glob import glob
import json
import pickle
from os.path import splitext
from os.path import basename	
from os.path import join


###################################

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/cobinding_analysis_refined"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

""" Choose the file category or TF category you want to analyse on """
file_category = "all_tf" #file_category = "all_tf"; file_category = "cr_tf"; file_category = "dbf_tf"

""" All TF motif containing dir """
input_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/fimo_motifs_for_hotanalysis"
dir_path = join(input_dir, "*motifs.bed")
fimo_motif_filelist = glob(dir_path) # represents all tf file list

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
for each_file in fimo_motif_filelist:
	tf_name = re.compile(r"(.*)_motifs.bed").findall(basename(each_file))[0]
	TF_check_list.append(tf_name)

# df_xls[df_xls["Target"].isin(TF_check_list)] # or,
if sorted(TF_check_list) == sorted(xls_tf_list):
	print "\nGreat!! TF list resembles b/w the file list and xls tf list..."
else:
	print "\nWarning: Your files TF list doesn't resemble the xls tf list..."


""" Select the category of TF for chosing the file list """
cr_tf_file_list = []
for each_tf in cr_df["Target"].tolist():
	for each_file in fimo_motif_filelist:
		if basename(each_file).startswith(each_tf):
			cr_tf_file_list.append(each_file)

cr_tf_file_list = list(set(cr_tf_file_list))

dbf_tf_file_list = []
for each_tf in dbf_df["Target"].tolist():
	for each_file in fimo_motif_filelist:
		if basename(each_file).startswith(each_tf):
			dbf_tf_file_list.append(each_file)

dbf_tf_file_list = list(set(dbf_tf_file_list))

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
		tf_name = re.compile(r"(.*)_motifs.bed").findall(basename(each_file))[0]
		df_read = pd.read_csv(each_file, header=None, sep="\t")
		df_read.columns = ["peak_chrom", "peak_start", "peak_end", "motif_id", "tf_name", "strand"]
		df_read["strand"] = "."
		df_read["tfmotif_id"] = df_read["tf_name"] + "|" + df_read["motif_id"]
		concat_file_list.append(df_read)
		print("Processing: {}\n".format(tf_name))

	combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
	sorted_combined_tf_df = combined_tf_df.sort_values(["peak_chrom", "peak_start", "peak_end"]).reset_index(drop=True)
	sorted_combined_tf_df.to_csv(join(output_dir, "final_combined_tf_coord_df.bed"), header=False, index=False, sep="\t")
	return(sorted_combined_tf_df)

# combined_tf_df = combine_tf_bedfiles(fimo_motif_filelist)


def cobinding_tf_peaks_count_perc_overlap_info(TF_file_list_of_interest, combined_TF_bedfile):
	# combined_TF_bedfile = combined_tf_df
	# TF_file_list_of_interest = dbf_tf_file_list
	# each_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/fimo_motifs_for_hotanalysis/JUN_human_motifs.bed"
	
	tf_file_list = TF_file_list_of_interest
	combined_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_TF_bedfile)
	master_tf_dict = {}
	master_check_list = []
	master_name = []

	for each_file in tf_file_list:
		tf_cobind_dict = {}
		tf_name = re.compile(r"(.*)_motifs.bed").findall(basename(each_file))[0]
		print "\nCurrently Processing %s\n" %(tf_name)
		df_read = pd.read_csv(each_file, sep="\t", header=None)
		df_read.columns = ["tf_chrom", "tf_start", "tf_end", "motif_id", "tf_name", "strand"]
		df_read["tfmotif_id"] = df_read["tf_name"] + "|" + df_read["motif_id"]
		# df_read["origin"] = df_read.apply(lambda row: "|".join(row.astype(str)), axis=1)
		
		pybed_df_list = []
		final_cobinding_tf_list = []
		cobinding_tf_peaks_count_list = []

		# Motifwise or Groupwise operation - treating each motif analysis cooccurrence as an independent event  
		for motif_id, df_motif in df_read.groupby(["motif_id"]):
			# df_motif = df_read.groupby(["motif_id"]).get_group("MOTIF4")
			df_motif_bedfile = pybedtools.BedTool.from_dataframe(df_motif)
			print df_motif_bedfile.head()

			"""Pybedtool intersection of tf_bedfile and combined_tf_bedfile"""
			print "Processing intersection for", tf_name
			tf_tf_intersect = df_motif_bedfile.intersect(combined_tf_bedfile, wa=True, wb=True)
			print tf_tf_intersect.head()

			pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
			header = df_motif.columns.tolist() + ["peak_chrom", "peak_start", "peak_end", "cobind_motif_id", "cobind_tf_name", "cobind_strand", "cobind_tfmotif_id"]
			pybed_df.columns = header

			"""Comment pybed_df if you need the original TF or the self_Tf-self_Tf intersect too"""
			pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] 
			cobinding_tf_list = pybed_df.loc[:,"cobind_tf_name"].unique().tolist()
			pybed_df_list.append(pybed_df)
			final_cobinding_tf_list.extend(cobinding_tf_list)
			
			# Add itself to the list of complex, but comment the line if you need self_TF-self_TF intersect
			# cobinding_tf_list.append(tf_name) # since; pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] commented;

			# Find unique non-repeating peaks wrt to chrom, start, tf_name, cobind_tf_name:  
			pybed_df_group = pybed_df.groupby(["tf_chrom","tf_start","tf_end","tfmotif_id","cobind_tfmotif_id", ]).count().reset_index()
			
			# debug = pybed_df_group[(pybed_df_group["tfmotif_id"] == "JUN_human|MOTIF4") & (pybed_df_group["cobind_tfmotif_id"] == "JUND_human|MOTIF2")]
			# debug = cobinding_tf_peaks_count[cobinding_tf_peaks_count["cobind_tfmotif_id"] == "JUND_human|MOTIF2"]

			cobinding_tf_peaks_count = pybed_df_group.groupby(["tfmotif_id","cobind_tfmotif_id"]).count().iloc[:,0].reset_index()		
			cobinding_tf_peaks_count = cobinding_tf_peaks_count.rename(columns={"tf_chrom" : "peak_overlap_count"})
			cobinding_tf_peaks_count["total_peak_count"] = df_motif.shape[0] # initial or original peak count
			cobinding_tf_peaks_count["percent_overlap"] = (cobinding_tf_peaks_count["peak_overlap_count"]/cobinding_tf_peaks_count["total_peak_count"])*100
			cobinding_tf_peaks_count = cobinding_tf_peaks_count.sort_values(["percent_overlap"], ascending=False)
			cobinding_tf_peaks_count_list.append(cobinding_tf_peaks_count)

		final_pybed_df = pd.concat(pybed_df_list, ignore_index=True)
		final_cobinding_tf_peaks_count = pd.concat(cobinding_tf_peaks_count_list, ignore_index=True)

		tf_cobind_dict["cobind_complex"] = final_pybed_df
		tf_cobind_dict["cobind_list"] = list(set(final_cobinding_tf_list))
		tf_cobind_dict["cobind_peak_count"] = final_cobinding_tf_peaks_count

		master_tf_dict[tf_name] = tf_cobind_dict
		master_check_list.append(cobinding_tf_list)
		master_name.append(tf_name)

	return(master_tf_dict)


""" Final Cobind peak count dataframe """
def cobinding_tf_peaks_overlap_final_output(master_tf_dict_info, File_category):
	master_TF_dict = master_tf_dict_info
	# master_TF_dict = master_TF_dict
	file_suffix =  File_category
	# file_suffix =  "dbf"
	cobind_peakcount_df_list = []
	for each_tf in master_TF_dict:
		cobind_peakcount_df_list.append(master_TF_dict[each_tf]["cobind_peak_count"])

	cobind_peakcount_dframe = pd.concat(cobind_peakcount_df_list)
	cobind_peak_count_df = cobind_peakcount_dframe.set_index(["tfmotif_id", "cobind_tfmotif_id"]) ### For heirarchial multiple indexing
	cobind_peak_count_df.columns = ["overlap_count", "total_peaks", "percent_overlap"]
	cobind_peak_count_df.to_csv(join(output_dir,"final_cobinding_peak_count_perc_overlap_" + file_suffix + ".bed"), sep="\t", header=True, index=True)

	cobind_peak_count_sorted = cobind_peak_count_df.sort_values(["percent_overlap"], ascending = False)
	cobind_peak_count_sorted.to_csv(join(output_dir,"final_cobinding_peak_count_perc_overlap_sorted_" + file_suffix + ".bed"), sep="\t", header=True, index=True)

	### Generating a heatmap data from the cobind peakcount dataframe:
	peak_table = cobind_peakcount_dframe.iloc[:,[0,1,4]]
	peak_heatmap_table = peak_table.pivot(index="tfmotif_id", columns="cobind_tfmotif_id", values = "percent_overlap")
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
"""For "test_tf" & "all_tf" file category, fimo_motif_filelist works or is common to both"""
combined_tf_df = combine_tf_bedfiles(dbf_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
combined_tf_df.to_pickle(join(output_dir, "final_combined_tf_coord_df.pkl"))
#combined_tf_df = combine_tf_bedfiles(cr_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
#combined_tf_df = combine_tf_bedfiles(fimo_motif_filelist) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist

master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(dbf_tf_file_list, combined_tf_df) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
#master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(cr_tf_file_list, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
#master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(fimo_motif_filelist, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist

final_cobind_df,final_cobind_df_sorted = cobinding_tf_peaks_overlap_final_output(master_TF_dict, file_category) 
final_cobind_df_sorted.to_pickle(join(output_dir, "final_cobind_df_sorted.pkl"))

# Merge and replace total motif peak count with motif specific peak count:
tfmotif_sp_count = combined_tf_df.groupby(["tfmotif_id"]).size().reset_index(name="motif_sp_count")
final_cobind_df_sorted_new = final_cobind_df_sorted.reset_index()

final_cobind_df_sorted_merged = pd.merge(final_cobind_df_sorted_new, tfmotif_sp_count,
											on = "tfmotif_id", how="left")
final_cobind_df_sorted_merged["percent_overlap"] = (final_cobind_df_sorted_merged["overlap_count"]/final_cobind_df_sorted_merged["motif_sp_count"])*100
final_cobind_df_sorted_merged_new = final_cobind_df_sorted_merged.sort_values(["percent_overlap"], ascending=False).reset_index(drop=True)

# To avoid confusion, dropping total peak count:
final_cobind_df_sorted_merged_new.drop(["total_peaks"], axis=1, inplace=True)
select_cols = ["tfmotif_id", "cobind_tfmotif_id", "overlap_count", "motif_sp_count", "percent_overlap"]
final_cobind_df_sorted_merged_new = final_cobind_df_sorted_merged_new.loc[:,select_cols]
final_cobind_df_sorted_merged_new.to_csv(join(output_dir,"final_cobinding_motif_perc_overlap_sorted_" + file_category + ".bed"), sep="\t", header=True, index=True)

print "Job completed successfully!!"




