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
# Coordinated for motif-cooccurrence analysis:
upstr_x1 = 1 
upstr_x2 = 40
downstr_x1 = 1
downstr_x2 = 40

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/cooccurrence_analysis"
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

#combined_tf_df = combine_tf_bedfiles_midpt_centered(dbf_tf_file_list)


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


def cooccurring_tf_peaks_perc_overlap_info(TF_file_list_of_interest, combined_TF_bedfile, upstr_x1, upstr_x2, downstr_x1, downstr_x2):
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
		df_read["origin"] = df_read.apply(lambda row: "|".join(row.astype(str)), axis=1)
		
		pybed_df_list = []
		final_cobinding_tf_list = []
		cobinding_tf_peaks_count_list = []

		# Motifwise or Groupwise operation - treating each motif analysis cooccurrence as an independent event  
		for motif_id, df_motif in df_read.groupby(["motif_id"]):
			# df_motif = df_read.groupby(["motif_id"]).get_group("MOTIF4")
			df_motif_upstr = df_motif.copy()
			df_motif_downstr = df_motif.copy()

			# Set the coordinates for 2 categories and concat `em:
			df_motif_upstr["tf_start"] = (df_motif["tf_start"] - upstr_x2) # upstr_x2 = 40
			df_motif_upstr["tf_end"] = (df_motif["tf_start"] - upstr_x1) + 1 # upstr_x1 =1, and + 1 for bedtool exclusivity
			df_motif_upstr["origin"] = df_motif["origin"]

			df_motif_downstr["tf_start"] = (df_motif["tf_end"] + downstr_x1) # downstr_x1 = 1
			df_motif_downstr["tf_end"] = (df_motif["tf_end"] + downstr_x2) + 1 # downstr_x2 = 40, and + 1 for bedtool exclusivity
			df_motif_downstr["origin"] = df_motif["origin"]

			df_motif_concat = pd.concat([df_motif_upstr, df_motif_downstr], ignore_index=True)
			df_motif_concat = df_motif_concat.sort_values(["tf_chrom", "tf_start", "tf_end"]).reset_index(drop=True)
			tf_concat_bedfile = pybedtools.BedTool.from_dataframe(df_motif_concat)
			print tf_concat_bedfile.head()

			"""Pybedtool intersection of tf_bedfile and combined_tf_bedfile"""
			print "Processing intersection for", tf_name
			orig_motif_bedfile = pybedtools.BedTool.from_dataframe(df_motif)
			nonoverlap_motif_pybed = combined_tf_bedfile.intersect(orig_motif_bedfile, v=True)
	
			tf_tf_intersect = tf_concat_bedfile.intersect(nonoverlap_motif_pybed, wa=True, wb=True)
			print tf_tf_intersect.head()

			pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
			header = df_motif.columns.tolist() + ["peak_chrom", "peak_start", "peak_end", "cobind_motif_id", "cobind_tf_name", "cobind_strand", "cobind_tfmotif_id"]
			# header = ["tf_chrom", "tf_start", "tf_end", "tf_name", "peak_chrom", "peak_start", "peak_end", "start", "end", "strand", "cobind_tf_name"]
			pybed_df.columns = header

			"""Comment pybed_df if you need the original TF or the self_Tf-self_Tf intersect too"""
			pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] 
			cobinding_tf_list = pybed_df.loc[:,"cobind_tf_name"].unique().tolist()
			pybed_df_list.append(pybed_df)
			final_cobinding_tf_list.extend(cobinding_tf_list)
			
			# Add itself to the list of complex, but comment the line if you need self_TF-self_TF intersect
			# cobinding_tf_list.append(tf_name) # since; pybed_df = pybed_df[~(pybed_df.loc[:,"cobind_tf_name"] == tf_name)] commented;

			# Find unique non-repeating peaks wrt to chrom, start, tf_name, cobind_tf_name:  
			pybed_df_group = pybed_df.groupby(["tf_chrom","tf_start","tf_end","tfmotif_id","origin","cobind_tfmotif_id", ]).count().reset_index()
			# debug = pybed_df_group[(pybed_df_group["tfmotif_id"] == "JUN_human|MOTIF4") & (pybed_df_group["cobind_tfmotif_id"] == "JUND_human|MOTIF2")]
			
			# Groupby tf_name, and cobind_tf_name to get total peak_overlap_count:
			# Find unique motif-origin - as motifs are broken up prior step; so we don't want same motif showing motif-occurrence with splitted motif-coord id:
			pybed_df_group_uniq = pybed_df_group.groupby(["origin","tfmotif_id","cobind_tfmotif_id"]).count().iloc[:,0].reset_index()
			cobinding_tf_peaks_count = pybed_df_group_uniq.groupby(["tfmotif_id","cobind_tfmotif_id"]).count().iloc[:,0].reset_index()		
			
			# debug = cobinding_tf_peaks_count[cobinding_tf_peaks_count["cobind_tfmotif_id"] == "JUND_human|MOTIF2"]
			cobinding_tf_peaks_count = cobinding_tf_peaks_count.rename(columns={"origin" : "peak_overlap_count"})
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
# file_category = "cr_tf" # cr_tf; dbf_tf, all_tf
# file_category = "all_tf" # cr_tf; dbf_tf, all_tf

output_dir = join(main_dir, file_category)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

"""Note: Make sure that tf_file_list matches with the file_category""" 
"""For "test_tf" & "all_tf" file category, fimo_motif_filelist works or is common to both"""
combined_tf_df = combine_tf_bedfiles(dbf_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
combined_tf_df.to_pickle(join(output_dir, "final_combined_tf_coord_df.pkl"))
# combined_tf_df = combine_tf_bedfiles(cr_tf_file_list) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
# combined_tf_df = combine_tf_bedfiles(fimo_motif_filelist) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist

master_TF_dict = cooccurring_tf_peaks_perc_overlap_info(dbf_tf_file_list, combined_tf_df, 1, 40, 1, 40) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
with open(join(output_dir, "master_TF_dict_motif_cooccurrence_final.pkl"), 'wb') as handle:
	pickle.dump(master_TF_dict, handle)
# master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(cr_tf_file_list, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist
# master_TF_dict = cobinding_tf_peaks_count_perc_overlap_info(fimo_motif_filelist, combined_tf_bedfile) # cr_tf_file_list; dbf_tf_file_list; fimo_motif_filelist

final_cobind_df, final_cobind_df_sorted = cobinding_tf_peaks_overlap_final_output(master_TF_dict, file_category) 
final_cobind_df_sorted.to_pickle(join(output_dir, "final_cobind_df_sorted.pkl"))

# Merge and replace total motif peak count with motif specific peak count:
tfmotif_sp_count = combined_tf_df.groupby(["tfmotif_id"]).size().reset_index(name="motif_sp_count")
final_cobind_df_sorted_new = final_cobind_df_sorted.reset_index()

final_cobind_df_sorted_merged =  pd.merge(final_cobind_df_sorted_new, tfmotif_sp_count,
											on = "tfmotif_id", how="left")
final_cobind_df_sorted_merged["percent_overlap"] = (final_cobind_df_sorted_merged["overlap_count"]/final_cobind_df_sorted_merged["motif_sp_count"])*100
final_cobind_df_sorted_merged_new = final_cobind_df_sorted_merged.sort_values(["percent_overlap"], ascending=False).reset_index(drop=True)

# To avoid confusion, dropping total peak count:
final_cobind_df_sorted_merged_new.drop(["total_peaks"], axis=1, inplace=True)
select_cols = ["tfmotif_id", "cobind_tfmotif_id", "overlap_count", "motif_sp_count", "percent_overlap"]
final_cobind_df_sorted_merged_new = final_cobind_df_sorted_merged_new.loc[:,select_cols]
final_cobind_df_sorted_merged_new.to_csv(join(output_dir,"final_cobinding_motif_perc_overlap_sorted_" + file_category + ".bed"), sep="\t", header=True, index=True)

print "Job completed successfully!!"


##############################################

##############################################

##############################################

import re, os 
from os.path import join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists

# Generate dict for peak-bed file list :
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"
bed_filelist = glob(file_pat)
bed_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_filedict = {bed_regex_pat.findall(file)[0]:file for file in bed_filelist}

# Generate dict for fimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_filelist = glob(file_pat)
fimo_regex_pat = re.compile(r'.*unique_TFs/.*_fimo_motif_(.*).txt')
fimo_filedict = {fimo_regex_pat.findall(file)[0]:file for file in fimo_filelist}


def final_fimo_motifs_model(fimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["motif_mid_point"] = (df["chromStart"] + df["chromEnd"])/2 
	df["motif_mid_point"] = df["motif_mid_point"].astype(int)
	df["motif_halfsize"] = (df["chromEnd"] - df["chromStart"])/2
	df["motif_halfsize"] = df["motif_halfsize"].astype(int)
    df["tf_name"] = tf_name 
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_mid_point", "motif_halfsize", "motif_id", "tf_name", "strand", "p-value"]
    motif_select_df = df.loc[:,select_cols]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)

    return(motif_sorted_df)


def final_fimo_motifs_model_with_pval(fimo_file, tf_name, pval, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["motif_mid_point"] = (df["chromStart"] + df["chromEnd"])/2 
	df["motif_mid_point"] = df["motif_mid_point"].astype(int)
	df["motif_halfsize"] = (df["chromEnd"] - df["chromStart"])/2
	df["motif_halfsize"] = df["motif_halfsize"].astype(int)
    df["tf_name"] = tf_name 
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_mid_point", "motif_halfsize", "motif_id", "tf_name", "strand", "p-value"]
    motif_select_df = df.loc[:,select_cols]
    #motif_select_df = motif_select_df.loc[motif_select_df["p-value"] < pval]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)

    return(motif_sorted_df)


def combine_tf_bedfiles_midpt_centered(TF_file, TF_name):
	tf_name = TF_name
	df_read = pd.read_csv(TF_file, header=None, sep="\t")
	df_read = df_read.iloc[:,[0,1,2]]
	df_read.columns = ["peak_chrom", "start", "end"]
	df_read["mid_point"] = (df_read["start"] + df_read["end"])/2 
	df_read["mid_point"] = df_read["mid_point"].astype(int)
	df_read["tf_name"] = tf_name
	select_cols = ["peak_chrom","peak_start","peak_end","start","end", "strand", "tf_name"]
	select_cols = ["peak_chrom","start","end", "mid_point", "tf_name"]
	df_read = df_read.loc[:,select_cols]
	df_final = df_read.sort_values(["peak_chrom","start","end"]).reset_index(drop=True)
	print df_final
	return(df_final)


############################
# Function implementations #
############################


output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

empty_fimo_coord_dict = {}
fimo_master_tfdict = {}
for tf in fimo_filedict:
    print("Currently processing : {} TF\n".format(TF_name))
    fimo_file = fimo_filedict["CTCF"]
    TF_name = "CTCF"
    tf = "CTCF"
    fimo_file = fimo_filedict[tf]
    TF_name = tf
    #fimo_coord_df = final_fimo_motifs_model_with_pval(fimo_file, TF_name, 1e-05, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)       
    fimo_coord_df = final_fimo_motifs_model(fimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)       

    if not fimo_coord_df.empty:
	    # Intersection of motif and total peak file:
	    motif_bedfile = pybedtools.BedTool.from_dataframe(fimo_coord_df)
	    peaks_df = combine_tf_bedfiles_midpt_centered(bed_filedict[TF_name], TF_name)
	    peaks_bedfile = pybedtools.BedTool.from_dataframe(peaks_df)
	    motif_peak_intersect = motif_bedfile.intersect(peaks_bedfile, wa=True, wb=True)
	    
	    # Read intersect file for further processing:
	    motif_peak_df = pd.read_csv(motif_peak_intersect.fn, sep="\t", header=None)
	    motif_peak_df.columns = fimo_coord_df.columns.tolist() + peaks_df.columns.tolist()
	    motif_peak_df["offset_dist"] = motif_peak_df["mid_point"] - motif_peak_df["motif_mid_point"]
		motif_peak_df["abs_offset_dist"] = abs(motif_peak_df["offset_dist"])

		# Find nearest motif to the peak
		
		gpd_df = motif_peak_df.groupby(["peak_chrom", "start", "end"])
		gpd_df_list = []
		for key,val in gpd_df:
			motif_peak_gp = gpd_df.get_group[key]
			motif_peak_gp = each_group.nsmallest(1, "abs_offset_dist")
			gpd_df_list.append(motif_peak_gp)
		motif_peak_gp_df = pd.concat(gpd_df_list, ignore_index=True)
		motif_offset_df = motif_peak_gp_df.groupby(["motif_id"])["abs_offset_dist"].median().reset_index()
		motif_size_df = motif_peak_gp_df.groupby(["motif_id"])["motif_halfsize"].median().reset_index()
		merged_motif_df = pd.merge(motif_offset_df, motif_size_df, on=["motif_id"])
		merged_motif_df["adj_offset_dist"] = merged_motif_df["abs_offset_dist"] - merged_motif_df["motif_halfsize"]
		merged_motif_df.loc[merged_motif_df["adj_offset_dist"] < 0, "adj_offset_dist"] = 0

		if TF_name not in fimo_master_tfdict:
			fimo_master_tfdict[TF_name] = merged_motif_df
	else:
		print("\n{} has no significant fimo coordinates".format(TF_name))
		empty_fimo_coord_dict[TF_name] = np.nan
		continue

# Concat all TF's data to compile merged fimo-tomtom-centrimo df:       
combined_motif_df = pd.concat(fimo_master_tfdict).reset_index().drop(["level_1"],axis=1)
combined_motif_df.rename(columns={"level_0":"TF_NAME", "motif_id" : "MOTIF_ID"},inplace=True)
combined_motif_df["TF_NAME_2"] = combined_motif_df["TF_NAME"].str.replace("_human", "").str.replace("\[FLAG\]", "").str.upper().str.replace(r"_.*", "")
combined_motif_df["MOTIF_ID_2"] = combined_motif_df["MOTIF_ID"].str.replace("MOTIF", "")
combined_motif_df["tf_motif_id"] = combined_motif_df["TF_NAME_2"] + "." + combined_motif_df["MOTIF_ID_2"]
combined_motif_df.to_pickle(join(output_dir,"combined_motif_offset_peak_df.pkl"))
combined_motif_df = pd.read_pickle(join(output_dir,"combined_motif_offset_peak_df.pkl"))


novel_filtered_motif_df = pd.read_csv(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50.txt"), sep="\t")
offset_motif_df = combined_motif_df[["tf_motif_id", "abs_offset_dist", "motif_halfsize", "adj_offset_dist"]]
final_merged_df1 = pd.merge(novel_filtered_motif_df, offset_motif_df, on=["tf_motif_id"])
final_merged_df.to_csv(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset.txt"), sep="\t", header=True, index=False)
