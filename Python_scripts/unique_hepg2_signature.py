import pandas as pd, numpy as np
import pybedtools
import scipy 
import os, re, glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from rpy2.robjects.packages import importr  #stats = importr('stats')
from rpy2.robjects.vectors import FloatVector
from os.path import join, basename, splitext
import pickle

# ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

""" Select the state number, to be analysed on """
#select_state_num = range(1,9+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
select_state_num = range(1,9) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:

""" Generate file list"""
ideas_file_dir = "/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state"

full_file_list = glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) # Includes HepG2 file
hepg2_file_list = [ each_file for each_file in glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) 
				if os.path.basename(each_file).startswith("HepG2")]
""" Excludes hepg2 file in the list"""
file_list = [ each_file for each_file in glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) 
						if not os.path.basename(each_file).startswith("HepG2")]

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/unique_hepg2_analysis_total"
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

	""" Read the mnemonics file and create a dictionary:
	# read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
	# state_list = read_file["Mnemonics"]
	# state_num_list = [ i for i in range(1,len(state_list)+1) ]
	# ideas_state_dict = dict(zip(state_num_list,state_list))
	# # print ideas_state_dict """

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
	#ideas_sorted_df.to_csv(join(output_dir, "final_ideas_prom_enh_concat_sorted.bed"))
	ideas_sorted_df.to_csv(join(output_dir, "final_ideas_prom_enh_concat_sorted_5celltype.bed"))

	""" Generate pybed object"""
	ideas_pybed = pybedtools.BedTool.from_dataframe(ideas_sorted_df) 
	merge_ideas_pybed = ideas_pybed.merge(c=[4,5], o=["distinct","distinct"])
	# merge_ideas_df = pd.read_csv(merge_ideas_pybed.fn, sep="\t")
	return(merge_ideas_pybed)

""" Merging the combined_ideas and hepg2 file by calling a function"""
hepg2_ideas_pybed =  merge_ideas_pybed(ideas_mnemonics_file, hepg2_file_list, select_state_num)
combined_ideas_pybed =  merge_ideas_pybed(ideas_mnemonics_file, file_list, select_state_num)
	

""" Final intersection b/w merged combined_ideas & merged hepg2_ideas"""
pybed_outfile_v = join(output_dir, "final_unique_hepg2_specific_regions_test.bed")
hepg2_unique_ideas_pybedfile = hepg2_ideas_pybed.intersect(combined_ideas_pybed, wa=True, wb=True, v = True, output=pybed_outfile_v)

""" Extract and count the number of different promoter and enhancer states """
uniq_sites_df = pd.read_csv(hepg2_unique_ideas_pybedfile.fn, sep="\t", header=None)
uniq_sites_df.columns = ["chrom", "start", "end", "state", "cellLine_state"]


""" Preparing for significant peak and background peaks overlap """
# dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)

null_generate_script = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_generate.py"
null_parameters = "-x 1 -r 1 "
null_hg19_indices = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_indices_hg19/"

master_ideas_dict = {}
master_sigcount_dict = {}
master_bgcount_dict = {}
total_peak_count_dict  = {}
total_uniqsite_count_dict = {}

for each_file in tf_file_list:
	tf_ideas_dict = {}
	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
	print "\nCurrently Processing %s\n" %(tf_name)
	df_read = pd.read_csv(each_file, sep="\t", header=None)
	df_read = df_read.iloc[:,[0,1,2]]
	df_read.columns = ["tf_chrom", "start", "end"]
	df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
	df_read["tf_start"] = df_read["mid_point"].astype(int) - (50)
	df_read["tf_end"] = df_read["mid_point"].astype(int) + (50)
	df_read["tf_name"] = tf_name
	total_peak_count_dict[tf_name] = df_read.shape[0]


	""" Preparing the file for significant peak overlap"""
	select_cols = ["tf_chrom","tf_start","tf_end","tf_name"]
	df_read = df_read.loc[:,select_cols]
	tf_bedfile = pybedtools.BedTool.from_dataframe(df_read)
	print "Preparing the intersection for", tf_name 
	print df_read


	""" Preparing the file for background peak overlap; bg stands for background"""
	#select_cols_bg = ["tf_chrom","tf_start","tf_end"]
	#df_bg_read = df_read.loc[:,select_cols_bg]
	#print df_bg_read
	df_read.to_csv(join(output_dir, "pre_background_file_"+tf_name), sep="\t", header=False, index=False)
	os.environ["null_generate_script"] = null_generate_script
	os.environ["null_parameters"] = null_parameters
	os.environ["null_hg19_indices"] = null_hg19_indices
	os.environ["input_file"] = join(output_dir, "pre_background_file_"+tf_name)
	os.environ["output_file"] = join(output_dir, "background_file_"+tf_name)
	
	if not os.path.exists(join(output_dir, "background_file_"+tf_name)):
		print "Generating the null sequence for", tf_name
		os.system("python $null_generate_script $null_parameters -o $output_file $input_file hg19 $null_hg19_indices")
	else:
		print "Background peak file exists for", tf_name
	#print "Preparing the background intersection for", tf_name 
	tf_bg_bedfile = pybedtools.BedTool(join(output_dir, "background_file_"+tf_name))
	print tf_bg_bedfile.head()

	
	""" Intersection b/w significant and each states, and likewise for backgound peaks """ 
	sig_count_list = []
	sig_state_list = []
	bg_count_list = []
	bg_state_list = []

	for each_state in target_state:
		each_uniq_df = uniq_sites_df[uniq_sites_df["state"].isin([each_state])]
		each_uniq_bedfile = pybedtools.BedTool.from_dataframe(each_uniq_df)
		print each_uniq_bedfile.head()
		if each_state not in total_uniqsite_count_dict:
			total_uniqsite_count_dict[each_state] = each_uniq_df.shape[0]

		""" For significant peak overlap"""
		print "\nCurrently processing the intersection for.... %s and %s" %(tf_name, each_state)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_state + "_unique_ideas_intersect.bed"))
		tf_ideas_intersect = tf_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_ideas_intersect.count() > 0:
			each_intersect_df = pd.read_csv(tf_ideas_intersect.fn, sep="\t", header=None)
			print tf_ideas_intersect.head()
		else:
			each_intersect_df = None

		### For the first time, need to create a dictionary, the update/merge the new dictionary into this
		### first time created dict: tf_ideas_dict[each_state] = {"sig_intersect": each_intersect_df}
		tf_ideas_dict[each_state] = {"sig_intersect": each_intersect_df}
		tf_ideas_dict[each_state].update({"sig_count": tf_ideas_intersect.count()})
		sig_count_list.append(tf_ideas_intersect.count())
		sig_state_list.append(each_state)


		""" For background peak overlap; bg stands for background"""		
		print "\nCurrently processing the background intersection for.... %s and %s" %(tf_name, each_state)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_state + "_background_intersect.bed"))
		tf_bg_ideas_intersect = tf_bg_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_bg_ideas_intersect.count() > 0:
			each_bg_intersect_df = pd.read_csv(tf_bg_ideas_intersect.fn, sep="\t", header=None)
			print tf_bg_ideas_intersect.head()
		else:
			each_bg_intersect_df = None
					
		tf_ideas_dict[each_state].update({"bg_intersect": each_bg_intersect_df})
		tf_ideas_dict[each_state].update({"bg_count": tf_bg_ideas_intersect.count()})
		bg_count_list.append(tf_bg_ideas_intersect.count())
		bg_state_list.append(each_state)

	master_sigcount_dict[tf_name] = dict(zip(sig_state_list, sig_count_list))
	master_bgcount_dict[tf_name] = dict(zip(bg_state_list, bg_count_list))
	master_ideas_dict[tf_name] = tf_ideas_dict

with open(join(output_dir,"final_master_ideas_dict.pkl"), "w") as outfile:
	pickle.dump(master_ideas_dict, outfile)

### To load back the pickled python objects:
# with open(join(output_dir,"final_master_ideas_dictas_dict_3tf.pkl")) as infile:
#	df = pickle.load(infile)

master_sig_bg_dict = {}
for each_tf in master_sigcount_dict:
	sig_key_list =  master_sigcount_dict[each_tf].keys()
	sig_value_list = master_sigcount_dict[each_tf].values()
	bg_key_list =  master_bgcount_dict[each_tf].keys()
	bg_value_list = master_bgcount_dict[each_tf].values()
	total_peak_count = total_peak_count_dict[each_tf]
	ideas_states_list = total_uniqsite_count_dict.keys()
	total_sites_list = total_uniqsite_count_dict.values()
	master_sig_bg_dict[each_tf]=pd.DataFrame({"sig_state":sig_key_list, "sig_hits":sig_value_list, "bg_state": bg_key_list, "bg_hits": bg_value_list, \
												 "total_peaks":total_peak_count, "states":ideas_states_list, "total_sites":total_sites_list})

combine_tf_df = pd.concat(master_sig_bg_dict)	
final_tf_df = combine_tf_df.reset_index().drop("level_1", axis=1)	
final_tf_df = final_tf_df.rename(columns={"level_0" : "tf_name"})
select_cols = ["tf_name", "states", "sig_hits", "bg_hits", "total_sites", "total_peaks"]
final_ideas_df = final_tf_df.loc[:,select_cols]
final_ideas_df.to_csv(join(output_dir,"final_tf_ideas_peak_count_with_bg.bed"), sep="\t", header=True, index=False)

with open(join(output_dir,"final_ideas_df.pkl"), "w") as outfile:
	pickle.dump(final_ideas_df, outfile)

with open(join(output_dir,"final_ideas_df.pkl")) as infile:
    final_ideas_df = pickle.load(infile)

final_ideas_df["sig_hits_unbound"] = final_ideas_df["total_peaks"] - final_ideas_df["sig_hits"]
final_ideas_df["bg_hits_unbound"] = final_ideas_df["total_peaks"] - final_ideas_df["bg_hits"]

# calculate fisher exact test & bh/bonferroni correction for multiple hypothesis/test correction:
def fishers_test(df_rows):
	x1 = df_rows["sig_hits"]
	x2 = df_rows["sig_hits_unbound"]
	y1 = df_rows["bg_hits"]
	y2 = df_rows["bg_hits_unbound"]
	pvalue = stats.fisher_exact([[x1,y1],[x2,y2]], "greater")[1] # R: 
	# pvalue = fisher.test(rbind(c(1,9),c(11,3)), alternative="greater")$p.value
	# by default the R's fisher exact has "lesser" as alt. hypothesis. So, to match
	# with python, set the alt hypothesis as; alternative = "less"

	""" #total_pval_test = df_rows["sig_hits"].shape[0]
		#bonferroni_correction = min(pvalue*total_pval_test, 1)
		#return(pvalue, df_rows["tf_name"]) """

	return pd.Series({"pvalue": pvalue, "sig_hits": x1, "bg_hits": y1, "tf_name": df_rows["tf_name"], 
						"states": df_rows["states"], "total_sites": df_rows["total_sites"], "total_peaks": df_rows["total_peaks"] })

final_ideas_pval_df = final_ideas_df.apply(fishers_test, axis=1)

# Replace the enrichment value 0 with the highest decimal place the python can generate i.e 1.0e-323
final_ideas_pval_df = final_ideas_pval_df.sort_values(["pvalue"])
final_ideas_pval_df["pvalue"] = final_ideas_pval_df["pvalue"].replace("0", 1.0e-323)

pvalue_list = final_ideas_pval_df["pvalue"].tolist()
rej, pval_corr = smm.multipletests(pvalue_list, alpha=0.05, method="fdr_bh")[:2] # method="bonferroni" or "hs"; if needed
final_ideas_pval_df["pval_corr"] = list(pval_corr)

order_cols = ["tf_name", "states", "sig_hits", "bg_hits", "pvalue", "pval_corr", "total_sites", "total_peaks"]
final_pval_qval_df =  final_ideas_pval_df.loc[:,order_cols]
final_pval_qval_df["percent_overlap"] = (final_pval_qval_df["sig_hits"]/final_pval_qval_df["total_peaks"])*100
pval_sorted_df = final_pval_qval_df.sort_values("pvalue")
sig_hits_sorted_df = final_pval_qval_df.sort_values("sig_hits", ascending=False)

# Log transformation:
final_pval_qval_df["-log10(pvalue)"] = -np.log10(final_pval_qval_df["pvalue"].astype(float))
final_pval_qval_df["-log10(pval_corr)"] = -np.log10(final_pval_qval_df["pval_corr"].astype(float))
final_pval_qval_df_sorted = final_pval_qval_df.replace(-0, 0)

with open(join(output_dir,"final_pval_qval_df.pkl"), "w") as outfile:
	pickle.dump(final_pval_qval_df_sorted, outfile)

with open(join(output_dir,"final_pval_qval_df.pkl")) as infile:
    final_pval_qval_df_sorted = pickle.load(infile)

""" Enrichment value of 0.001 """
final_pval_qval_001 = final_pval_qval_df_sorted[final_pval_qval_df_sorted["pval_corr"] < 0.001]
final_pval_qval_001.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.001.bed"), sep="\t", header=True, index=False)
final_pval_qval_heatmap_001 = final_pval_qval_001.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
final_pval_qval_heatmap_001.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.001.bed"), sep="\t", header=True, index=True)

# Based on the output of heatmap decided to exclude some of the datas:
final_pval_qval_state_filtered_001 = final_pval_qval_001[final_pval_qval_001["states"].isin(["Enh", "EnhF", "PromCtcf"])]
final_pval_qval_state_filtered_001.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.001_state_filtered.bed"), sep="\t", header=True, index=True)
final_pval_qval_heatmap_001 = final_pval_qval_state_filtered_001.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
final_pval_qval_heatmap_001.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.001_state_filtered.bed"), sep="\t", header=True, index=True)

""" Enrichmnent value of 0.01 """
final_pval_qval_01 = final_pval_qval_df_sorted[final_pval_qval_df_sorted["pval_corr"] < 0.01]
final_pval_qval_01.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.01.bed"), sep="\t", header=True, index=False)
final_pval_qval_heatmap_01 = final_pval_qval_01.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
final_pval_qval_heatmap_01.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.01.bed"), sep="\t", header=True, index=True)

# Based on the output of heatmap decided to exclude some of the datas:
final_pval_qval_state_filtered_01 = final_pval_qval_01[final_pval_qval_01["states"].isin(["Enh", "EnhF", "PromCtcf"])]
final_pval_qval_state_filtered_01.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.01_state_filtered.bed"), sep="\t", header=True, index=True)
final_pval_qval_heatmap_01 = final_pval_qval_state_filtered_01.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
final_pval_qval_heatmap_01.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.01_state_filtered.bed"), sep="\t", header=True, index=True)

### Herebelow codes are not needed for this particular analysis.

""" All TF peaks/binding sites files """
#dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)

master_tf_dict = {}

for each_file in tf_file_list:
	tf_ideas_dict = {}
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
	
	""" Pybedtool intersection of tf_bedfile and hepg2_unique_specific_bedfile """
	print "Processing the intersection for", tf_name
	pybed_outfile = join(output_dir, (tf_name + "_and_unique_hepg2_ideas_intersect.bed"))
	tf_ideas_intersect = tf_bedfile.intersect(hepg2_unique_ideas_pybedfile, wa=True, wb=True, output=pybed_outfile)
	print tf_ideas_intersect.head()

	if tf_ideas_intersect.count() > 0:
		pybed_df = pd.read_csv(tf_ideas_intersect.fn, sep="\t", header=None)
		header = ["tf_chrom", "tf_start", "tf_end", "tf_name", "chrom", "start", "end", "state", "state_cellLine"]
		pybed_df.columns = header

		# Find unique non-repeating peaks wrt to chrom, start, tf_name, cobind_tf_name:  
		pybed_df_group = pybed_df.groupby(["tf_name","state"]).count().reset_index()
		print pybed_df_group

		# Groupby tf_name, and cobind_tf_name to get total peak_overlap_count:
		tf_ideas_peaks_count = pybed_df_group.iloc[:,[0,1,2]]
		tf_ideas_peaks_count = tf_ideas_peaks_count.rename(columns={"tf_chrom" : "peak_overlap_count"})
		tf_ideas_peaks_count["total_peak_count"] = df_read.shape[0] # initial or original peak count
		tf_ideas_peaks_count["percent_overlap"] = (tf_ideas_peaks_count["peak_overlap_count"]/tf_ideas_peaks_count["total_peak_count"])*100
		tf_ideas_peaks_count = tf_ideas_peaks_count.sort_values(["percent_overlap"], ascending=False)

		tf_ideas_dict["tf_ideas_intersect"] = pybed_df
		tf_ideas_dict["tf_ideas_peak_count"] = tf_ideas_peaks_count

		if tf_name not in master_tf_dict:
			master_tf_dict[tf_name] = tf_ideas_dict


""" Final tf-ideas_state peak count dataframe on unique hepg2 """
tf_ideas_df_list = []
for each_tf in master_tf_dict:
	tf_ideas_df_list.append(master_tf_dict[each_tf]["tf_ideas_peak_count"])

tf_ideas_peakcount_dframe = pd.concat(tf_ideas_df_list)


tf_ideas_peak_count_df = tf_ideas_peakcount_dframe.set_index(["tf_name", "state"]) ### For heirarchial multiple indexing
tf_ideas_peak_count_df.columns = ["overlap_count", "total_peaks", "percent_overlap"]
tf_ideas_peak_count_df.to_csv(join(output_dir,"final_tf_ideas_peak_count.bed"), sep="\t", header=True, index=True)

tf_ideas_peak_count_sorted = tf_ideas_peak_count_df.sort_values(["percent_overlap"], ascending = False)
tf_ideas_peak_count_sorted.to_csv(join(output_dir,"final_tf_ideas_peak_count_sorted.bed"), sep="\t", header=True, index=True)

""" Generating a heatmap data from tf-ideas_state peak count dataframe """
peak_table = tf_ideas_peakcount_dframe.iloc[:,[0,1,4]]
peak_heatmap_table = peak_table.pivot(index="tf_name", columns="state", values = "percent_overlap")
peak_heatmap_table.to_csv(join(output_dir,"final_tf_ideas_peak_count_heatmap.bed"), sep="\t", header=True, index=True)

