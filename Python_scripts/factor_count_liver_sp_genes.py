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


main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_factor_count_with_ideas_bins"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

sp_gene_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/hepg2_liver_gene_sp_analysis/files_2kb_from_tss"

""" Read liver and hepg2 specific genes """
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Liver_specific_genes_HepG2.xlsx"
df_xls = pd.ExcelFile(excelfile).parse()
sp_gene_list = df_xls["Gene"].unique().tolist()
df_xls.columns = [ str(each).replace(" ", "_") for each in df_xls.columns.tolist()]
df_xls.columns
tpm_df = df_xls.loc[:, ["Gene", "HepG2_RNAseq", "Liver_RNAseq"]]
tpm_df.columns = ["gene_name", "HepG2_RNAseq", "Liver_RNAseq"]
	
sp_gene_df =  pd.read_csv(join(sp_gene_file_dir, "final_tf_gene_enrichment_with_pval_qval.bed"), sep="\t")
sp_gene_df = sp_gene_df.iloc[:, [0,1,2]]
sp_gene_df["bool_value"] = np.where(sp_gene_df["sig_hits"] > 0, 1, 0)
sp_gene_df.to_csv(join(sp_gene_file_dir, "final_tf_binding_with_gene_bins_boolean_style.bed"), sep="\t", header=True, index=False)
sp_gene_df.to_csv(join(main_dir, "final_tf_binding_with_gene_bins_boolean_style.bed"), sep="\t", header=True, index=False)

""" Merge the genes with tf count and TPM value to generate a barplot data """
tf_count_with_genes_bin = sp_gene_df.groupby(["gene_name"]).apply(lambda x : x["bool_value"].sum())
tf_count_with_genes_bin = tf_count_with_genes_bin.reset_index().rename(columns={0: "tf_counts"})
merged_gene_df = pd.merge(tf_count_with_genes_bin, tpm_df, left_on = "gene_name", right_on= "gene_name")
merged_gene_df = merged_gene_df.sort_values(["tf_counts","HepG2_RNAseq"])
merged_gene_df.to_csv(join(sp_gene_file_dir, "final_tf_binding_with_gene_bins_merged_barplot_data.bed"), sep="\t", header=True, index=False)
merged_gene_df.to_csv(join(main_dir, "final_tf_binding_with_gene_bins_merged_barplot_data.bed"), sep="\t", header=True, index=False)

""" Generate a heatmap data for the TF binding in each bin"""
sp_gene_df_heatmap_data = sp_gene_df.pivot(index="tf_name", columns="gene_name", values= "bool_value")
sp_gene_df_heatmap_data.to_csv(join(sp_gene_file_dir, "final_tf_binding_with_gene_bins_heatmap_data.bed"), sep="\t", header=True, index=True)
sp_gene_df_heatmap_data.to_csv(join(main_dir, "final_tf_binding_with_gene_bins_heatmap_data.bed"), sep="\t", header=True, index=True)


""" Extract the tf binding data for diff gene including ALB, APOA2 """
gene_df =  pd.read_csv(join(sp_gene_file_dir, "final_tf_gene_enrichment_with_pval_qval.bed"), sep="\t")
alb_tf_df = gene_df[(gene_df["gene_name"] == "ALB") & (gene_df["sig_hits"] >0)]
final_alb_tf_df = alb_tf_df.loc[:,["tf_name", "gene_name", "sig_hits"]] 
final_alb_tf_df.to_csv(join(sp_gene_file_dir, "ALB_promoter_tf_binding_data.txt"), sep="\t", header=True, index=False)

apoa2_tf_df = gene_df[(gene_df["gene_name"] == "APOA2") & (gene_df["sig_hits"] >0)] 
final_apoa2_tf_df = apoa2_tf_df.loc[:,["tf_name", "gene_name", "sig_hits"]] 
final_apoa2_tf_df.to_csv(join(sp_gene_file_dir, "APOA2_promoter_tf_binding_data.txt"), sep="\t", header=True, index=False)
