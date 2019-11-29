import pandas as pd
import numpy as np
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

null_sequence_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/unique_hepg2_analysis_total"
main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/hepg2_liver_gene_sp_analysis"
# main_dir = os.path.expanduser("~/Dropbox/for_chris/batch_I/hepg2_liver_gene_sp_analysis")
output_dir =  join(main_dir, "files_potential_enhancer")

if not os.path.exists(main_dir):
	os.makedirs(main_dir)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

refgene_file = "/gpfs/gpfs1/home/schhetri/for_chris/refGene_hg19"
#refgene_file = os.path.expanduser("~/Dropbox/for_chris/refGene_hg19")

""" All TF containing dir """
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
# dir_path = os.path.expanduser("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/idr_passed_peaks_total/unique_TFs/SL*")
all_tf_file_list = glob.glob(dir_path) # represents all tf file list


""" Read liver and hepg2 specific genes """
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Liver_specific_genes_in_HepG2_with_enh_and_tissue_specific_region.xlsx" 
# excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Liver_specific_genes_HepG2.xlsx"
# excelfile = os.path.expanduser("~/Dropbox/for_chris/Liver_specific_genes_HepG2.xlsx")
df_xls = pd.ExcelFile(excelfile).parse()
df_xls =  df_xls.iloc[:,[0,1]]
df_xls.columns = ["gene_name", "potential_enh_prom"]
# df_xls["splitted"] = df_xls["potential_enh_prom"].apply(lambda x: re.split("[:-]", x))
df_xls["splitted"] = [re.split("[: -]",each) for each in df_xls["potential_enh_prom"]]
df_xls[["chrom", "start", "end"]] = df_xls["splitted"].apply(pd.Series)
select_cols = ["chrom", "start", "end", "gene_name"]
df_xls = df_xls.loc[:,select_cols]

# some wierd float thing happened, so replaced it by the correct int format:
df_xls["end"].replace("94822533.", "94822533", inplace=True)
enh_coord_df = df_xls.sort_values(["chrom","start"])
final_sp_gene_list = enh_coord_df["gene_name"].unique().tolist()


""" Preparing for significant peak and background peaks overlap """
# dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)

null_generate_script = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_generate.py"
null_parameters = "-x 1 -r 1 "
null_hg19_indices = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_indices_hg19/"

master_gene_dict = {}
master_sigcount_dict = {}
master_bgcount_dict = {}
total_peak_count_dict  = {}
total_gene_dict = {}


for each_file in tf_file_list:
	tf_gene_dict = {}
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
	df_read.to_csv(join(null_sequence_dir, "pre_background_file_" + tf_name), sep="\t", header=False, index=False)
	os.environ["null_generate_script"] = null_generate_script
	os.environ["null_parameters"] = null_parameters
	os.environ["null_hg19_indices"] = null_hg19_indices
	os.environ["input_file"] = join(null_sequence_dir, "pre_background_file_"+tf_name)
	os.environ["output_file"] = join(null_sequence_dir, "background_file_"+tf_name)
	

	if not os.path.exists(join(null_sequence_dir, "background_file_"+tf_name)):
		print "Generating the null sequence for", tf_name
		os.system("python $null_generate_script $null_parameters -o $output_file $input_file hg19 $null_hg19_indices")
	else:
		print "Background peak file exists for", tf_name
	#print "Preparing the background intersection for", tf_name 
	tf_bg_bedfile = pybedtools.BedTool(join(null_sequence_dir, "background_file_"+tf_name))
	print tf_bg_bedfile.head()

	
	""" Intersection b/w significant and each states, and likewise for backgound peaks """ 
	sig_count_list = []
	sig_state_list = []
	bg_count_list = []
	bg_state_list = []

	for each_gene in final_sp_gene_list:
		each_uniq_df = enh_coord_df[enh_coord_df["gene_name"].isin([each_gene])]
		each_uniq_bedfile = pybedtools.BedTool.from_dataframe(each_uniq_df)
		print each_uniq_bedfile.head()

		if each_gene not in total_gene_dict:
			total_gene_dict[each_gene] = each_uniq_df.shape[0]

		""" For significant peak overlap"""
		print "\nCurrently processing the intersection for.... %s and %s" %(tf_name, each_gene)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_gene + "_unique_ideas_intersect.bed"))
		tf_gene_intersect = tf_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_gene_intersect.count() > 0:
			each_intersect_df = pd.read_csv(tf_gene_intersect.fn, sep="\t", header=None)
			print tf_gene_intersect.head()
		else:
			each_intersect_df = None

		### For the first time, need to create a dictionary, the update/merge the new dictionary into this
		### first time created dict: tf_gene_dict[each_gene] = {"sig_intersect": each_intersect_df}
		tf_gene_dict[each_gene] = {"sig_intersect": each_intersect_df}
		tf_gene_dict[each_gene].update({"sig_count": tf_gene_intersect.count()})
		sig_count_list.append(tf_gene_intersect.count())
		sig_state_list.append(each_gene)


		""" For background peak overlap; bg stands for background"""		
		print "\nCurrently processing the background intersection for.... %s and %s" %(tf_name, each_gene)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_gene + "_background_intersect.bed"))
		tf_bg_gene_intersect = tf_bg_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_bg_gene_intersect.count() > 0:
			each_bg_intersect_df = pd.read_csv(tf_bg_gene_intersect.fn, sep="\t", header=None)
			print tf_bg_gene_intersect.head()
		else:
			each_bg_intersect_df = None
					
		tf_gene_dict[each_gene].update({"bg_intersect": each_bg_intersect_df})
		tf_gene_dict[each_gene].update({"bg_count": tf_bg_gene_intersect.count()})
		bg_count_list.append(tf_bg_gene_intersect.count())
		bg_state_list.append(each_gene)

	master_sigcount_dict[tf_name] = dict(zip(sig_state_list, sig_count_list))
	master_bgcount_dict[tf_name] = dict(zip(bg_state_list, bg_count_list))
	master_gene_dict[tf_name] = tf_gene_dict

with open(join(output_dir,"final_master_gene_dict.pkl"), "w") as outfile:
	pickle.dump(master_gene_dict, outfile)

### To load back the pickled python objects:
# with open(join(output_dir,"final_master_gene_dictas_dict_3tf.pkl")) as infile:
#	df = pickle.load(infile)
master_sig_bg_dict = {}
for each_tf in master_sigcount_dict:
	sig_key_list =  master_sigcount_dict[each_tf].keys()
	sig_value_list = master_sigcount_dict[each_tf].values()
	bg_key_list =  master_bgcount_dict[each_tf].keys()
	bg_value_list = master_bgcount_dict[each_tf].values()
	total_peak_count = total_peak_count_dict[each_tf]
	total_gene_list = total_gene_dict.keys()
	total_sites_list = total_gene_dict.values()
	master_sig_bg_dict[each_tf]=pd.DataFrame({"sig_state":sig_key_list, "sig_hits":sig_value_list, "bg_state": bg_key_list, "bg_hits": bg_value_list, \
												 "total_peaks":total_peak_count, "gene_name":total_gene_list, "total_gene_sites": total_gene_dict.values()})

combine_tf_df = pd.concat(master_sig_bg_dict)	
final_tf_df = combine_tf_df.reset_index().drop("level_1", axis=1)	
final_tf_df = final_tf_df.rename(columns={"level_0" : "tf_name"})
select_cols = ["tf_name", "gene_name", "sig_hits", "bg_hits", "total_gene_sites", "total_peaks"]
final_ideas_df = final_tf_df.loc[:,select_cols]
final_ideas_df.to_csv(join(output_dir,"final_tf_gene_enh_peak_count_with_bg.bed"), sep="\t", header=True, index=False)

with open(join(output_dir,"final_sp_gene_df.pkl"), "w") as outfile:
	pickle.dump(final_ideas_df, outfile)

with open(join(output_dir,"final_sp_gene_df.pkl")) as infile:
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
						"gene_name": df_rows["gene_name"], "total_gene_sites": df_rows["total_gene_sites"], "total_peaks": df_rows["total_peaks"] })

final_ideas_pval_df = final_ideas_df.apply(fishers_test, axis=1)

# Replace the enrichment value 0 with the highest decimal place the python can generate i.e 1.0e-323
final_ideas_pval_df = final_ideas_pval_df.sort_values(["pvalue"])
final_ideas_pval_df["pvalue"] = final_ideas_pval_df["pvalue"].replace("0", 1.0e-323)

pvalue_list = final_ideas_pval_df["pvalue"].tolist()
rej, pval_corr = smm.multipletests(pvalue_list, alpha=0.05, method="fdr_bh")[:2] # method="bonferroni" or "hs"; if needed
final_ideas_pval_df["pval_corr"] = list(pval_corr)

order_cols = ["tf_name", "gene_name", "sig_hits", "bg_hits", "pvalue", "pval_corr", "total_gene_sites", "total_peaks"]
final_pval_qval_df =  final_ideas_pval_df.loc[:,order_cols]

# Log transformation:
final_pval_qval_df["-log2(pvalue)"] = -np.log2(final_pval_qval_df["pvalue"].astype(float))
final_pval_qval_df["-log2(pval_corr)"] = -np.log2(final_pval_qval_df["pval_corr"].astype(float))
final_pval_qval_df_sorted = final_pval_qval_df.replace(-0, 0)

with open(join(output_dir,"final_pval_qval_df.pkl"), "w") as outfile:
	pickle.dump(final_pval_qval_df_sorted, outfile)

with open(join(output_dir,"final_pval_qval_df.pkl")) as infile:
    final_pval_qval_df_sorted = pickle.load(infile)


""" Enrichmnent value heatmap """
final_pval_qval_df_sorted.to_csv(join(output_dir,"final_tf_gene_enh_enrichment_with_pval_qval.bed"), sep="\t", header=True, index=False)
final_pval_qval_heatmap_01 = final_pval_qval_df_sorted.pivot(index="tf_name", columns="gene_name", values = "-log2(pvalue)")
final_pval_qval_heatmap_01.to_csv(join(output_dir,"final_tf_gene_enrichment_pval_corr_heatmap_data.bed"), sep="\t", header=True, index=True)

# Based on the output of heatmap decided to exclude some of the datas:
sig_gene_list = final_pval_qval_df[final_pval_qval_df["pvalue"] < 0.05]["gene_name"].unique().tolist()
final_pval_qval_gene_filtered = final_pval_qval_df_sorted[final_pval_qval_df_sorted["gene_name"].isin(sig_gene_list)]
final_pval_qval_gene_filtered.to_csv(join(output_dir,"final_tf_gene_enh_enrichment_with_pval_qval_gene_filtered.bed"), sep="\t", header=True, index=True)
final_pval_qval_heatmap = final_pval_qval_gene_filtered.pivot(index="tf_name", columns="gene_name", values = "-log2(pvalue)")
final_pval_qval_heatmap.to_csv(join(output_dir,"final_tf_gene_enh_enrichment_pval_corr_heatmap_data_gene_filtered.bed"), sep="\t", header=True, index=True)

