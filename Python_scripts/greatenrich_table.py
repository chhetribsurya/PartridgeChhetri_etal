import pandas as pd
from scipy import stats
import os
import re
from os.path import join
from os.path import basename

output_dir = os.path.expanduser("~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/great_analysis/downloaded_final_data")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

""" Data from genomic_coord <- gene association table from great analysis for TF bound sites """
df = pd.read_csv(join(output_dir, "TF_bound_sites_promoter_gene_association.txt"), sep="\t", header=None)
df_prom = df.loc[~(df.iloc[:,1] == "NONE")]
df_prom = df_prom.iloc[1:,[0,1]] # gettting rid of 0th row(with no info)
df_prom.columns = ["genomic_coords", "gene_names"]
#df_final = pd.DataFrame(df_prom["gene_names"].str.split(',').tolist(), index=df_prom["genomic_coords"]).stack().reset_index().rename(columns={0:"gene_names"})
#df_final3 = df_prom.set_index(["genomic_coords"])["gene_names"].str.split(",",expand=True).stack().reset_index().rename(columns={0:'gene_names'}).loc[:, df_prom.columns]

""" Provide the name of column to be splitted by delimiter ",": (replace "gene_names" w the variable to be splitted) 
    Especially useful for the refseq gene splitting of exons and introns: (https://stackoverflow.com/questions/12680754/split-pandas-dataframe-string-entry-to-separate-rows)"""
df_prom_final = (df_prom.set_index(df_prom.columns.drop('gene_names',1).tolist())["gene_names"].str.split(',', expand=True).stack().reset_index().rename(columns={0:'gene_names'}).loc[:, df_prom.columns])

"""Further splitting of the column wrt delimiter "space" i.e default"""
df_prom_final_tf_bound = df_prom_final.set_index(["genomic_coords"])["gene_names"].str.split(expand=True).reset_index()
df_prom_final_tf_bound.columns = ["genomic_coords", "gene_names", "genomic_dist"]

"""Write the file"""
df_prom_final_tf_bound.to_csv(join(output_dir, "TF_bound_sites_promoter_gene_association_w_genomic_distance.txt"), sep="\t", header=True, index=True)


#####################################

""" Read the rna seq cell line data from protein atlas for HepG2 """
df = pd.read_csv(join(output_dir, "rna_celline.tsv"), sep="\t")
df_rnaseq = df.loc[df["Sample"] == "Hep G2"]

#####################################

""" Now merging the genes list associated to TF bound and unbound for promoter associated regions with the gene list from the HepG2 rnaseq to get the TPM values for each genes in common b/w them"""
df_prom_rnaseq_gene_merged_tf_bound = pd.merge(df_prom_final_tf_bound, df_rnaseq, left_on = "gene_names", right_on ="Gene name")
df_prom_rnaseq_gene_merged_tf_bound_0tpm = df_prom_rnaseq_gene_merged_tf_bound.loc[df_prom_rnaseq_gene_merged_tf_bound["Value"]  ==  0]
zero_tpm_perc =  df_prom_rnaseq_gene_merged_tf_bound_0tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_bound.shape[0])*100; zero_tpm_perc
df_prom_rnaseq_gene_merged_tf_bound_0tpm.shape[0], df_prom_rnaseq_gene_merged_tf_bound.shape[0], zero_tpm_perc

df_prom_rnaseq_gene_merged_tf_bound_lessthan_5tpm = df_prom_rnaseq_gene_merged_tf_bound.loc[df_prom_rnaseq_gene_merged_tf_bound["Value"]  < 5]
df_prom_rnaseq_gene_merged_tf_bound_lessthan_5tpm.shape[0]
lt_five_tpm_perc =  df_prom_rnaseq_gene_merged_tf_bound_lessthan_5tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_bound.shape[0])*100; lt_five_tpm_perc
df_prom_rnaseq_gene_merged_tf_bound_lessthan_5tpm.shape[0], df_prom_rnaseq_gene_merged_tf_bound.shape[0], lt_five_tpm_perc

df_prom_rnaseq_gene_merged_tf_bound_morethan_5tpm = df_prom_rnaseq_gene_merged_tf_bound.loc[df_prom_rnaseq_gene_merged_tf_bound["Value"]  >= 5]
egt_five_tpm_perc =  df_prom_rnaseq_gene_merged_tf_bound_morethan_5tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_bound.shape[0])*100; egt_five_tpm_perc
df_prom_rnaseq_gene_merged_tf_bound_morethan_5tpm.shape[0], df_prom_rnaseq_gene_merged_tf_bound.shape[0], egt_five_tpm_perc

""" Write the files """
df_prom_rnaseq_gene_merged_tf_bound.to_csv(join(output_dir, "TF_bound_prom_gene_assoc_patlas_hepg2_rna_seq_tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_bound_0tpm.to_csv(join(output_dir, "TF_bound_prom_gene_assoc_patlas_hepg2_rna_seq_0tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_bound_lessthan_5tpm.to_csv(join(output_dir, "TF_bound_prom_gene_assoc_patlas_hepg2_rna_seq_lessthan_5tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_bound_morethan_5tpm.to_csv(join(output_dir, "TF_bound_prom_gene_assoc_patlas_hepg2_rna_seq_equal_or_morethan_5tpm_merged.txt"), sep="\t", header=True, index=True)

############################################

""" Data from genomic_coord <- gene association table from great analysis for TF unbound sites """
df = pd.read_csv(join(output_dir, "TF_unbound_sites_gene_association.txt"), sep="\t", header=None)
df_prom = df.loc[~(df.iloc[:,1] == "NONE")]
df_prom = df_prom.iloc[1:,[0,1]] # gettting rid of 0th row(with no info)
df_prom.columns = ["genomic_coords", "gene_names"]
#df_final = pd.DataFrame(df_prom["gene_names"].str.split(',').tolist(), index=df_prom["genomic_coords"]).stack().reset_index().rename(columns={0:"gene_names"})
#df_final3 = df_prom.set_index(["genomic_coords"])["gene_names"].str.split(",",expand=True).stack().reset_index().rename(columns={0:'gene_names'}).loc[:, df_prom.columns]

""" Provide the name of column to be splitted by delimiter ",": (replace "gene_names" w the variable to be splitted) 
    Especially useful for the refseq gene splitting of exons and introns: (https://stackoverflow.com/questions/12680754/split-pandas-dataframe-string-entry-to-separate-rows) """
df_prom_final = (df_prom.set_index(df_prom.columns.drop('gene_names',1).tolist())["gene_names"].str.split(',', expand=True).stack().reset_index().rename(columns={0:'gene_names'}).loc[:, df_prom.columns])

""" Further splitting of the column wrt delimiter "space" i.e default """
df_prom_final_tf_unbound = df_prom_final.set_index(["genomic_coords"])["gene_names"].str.split(expand=True).reset_index()
df_prom_final_tf_unbound.columns = ["genomic_coords", "gene_names", "genomic_dist"]

""" Write the file """	
df_prom_final_tf_unbound.to_csv(join(output_dir, "TF_unbound_sites_promoter_gene_association_w_genomic_distance_1.txt"), sep="\t", header=True, index=True)


#####################################

""" Read the rna seq cell line data from protein atlas for HepG2 """
df = pd.read_csv(join(output_dir, "rna_celline.tsv"), sep="\t")
df_rnaseq = df.loc[df["Sample"] == "Hep G2"]

#####################################

""" Now merging the genes list associated to TF bound and unbound for promoter associated regions with the gene list from the HepG2 rnaseq to get the TPM values for each genes in common b/w them"""
df_prom_rnaseq_gene_merged_tf_unbound = pd.merge(df_prom_final_tf_unbound, df_rnaseq, left_on = "gene_names", right_on ="Gene name")
df_prom_rnaseq_gene_merged_tf_unbound_0tpm = df_prom_rnaseq_gene_merged_tf_unbound.loc[df_prom_rnaseq_gene_merged_tf_unbound["Value"]  ==  0]
zero_tpm_perc =  df_prom_rnaseq_gene_merged_tf_unbound_0tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_unbound.shape[0])*100; zero_tpm_perc
df_prom_rnaseq_gene_merged_tf_unbound_0tpm.shape[0], df_prom_rnaseq_gene_merged_tf_unbound.shape[0], zero_tpm_perc

df_prom_rnaseq_gene_merged_tf_unbound_lessthan_5tpm = df_prom_rnaseq_gene_merged_tf_unbound.loc[df_prom_rnaseq_gene_merged_tf_unbound["Value"]  < 5]
df_prom_rnaseq_gene_merged_tf_unbound_lessthan_5tpm.shape[0]
lt_five_tpm_perc =  df_prom_rnaseq_gene_merged_tf_unbound_lessthan_5tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_unbound.shape[0])*100; lt_five_tpm_perc
df_prom_rnaseq_gene_merged_tf_unbound_lessthan_5tpm.shape[0], df_prom_rnaseq_gene_merged_tf_unbound.shape[0], lt_five_tpm_perc

df_prom_rnaseq_gene_merged_tf_unbound_morethan_5tpm = df_prom_rnaseq_gene_merged_tf_unbound.loc[df_prom_rnaseq_gene_merged_tf_unbound["Value"]  >= 5]
egt_five_tpm_perc =  df_prom_rnaseq_gene_merged_tf_unbound_morethan_5tpm.shape[0]/float(df_prom_rnaseq_gene_merged_tf_unbound.shape[0])*100; egt_five_tpm_perc
df_prom_rnaseq_gene_merged_tf_unbound_morethan_5tpm.shape[0], df_prom_rnaseq_gene_merged_tf_unbound.shape[0], egt_five_tpm_perc

""" Write the files """
df_prom_rnaseq_gene_merged_tf_unbound.to_csv(join(output_dir, "TF_unbound_prom_gene_assoc_patlas_hepg2_rna_seq_tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_unbound_0tpm.to_csv(join(output_dir, "TF_unbound_prom_gene_assoc_patlas_hepg2_rna_seq_0tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_unbound_lessthan_5tpm.to_csv(join(output_dir, "TF_unbound_prom_gene_assoc_patlas_hepg2_rna_seq_lessthan_5tpm_merged.txt"), sep="\t", header=True, index=True)
df_prom_rnaseq_gene_merged_tf_unbound_morethan_5tpm.to_csv(join(output_dir, "TF_unbound_prom_gene_assoc_patlas_hepg2_rna_seq_equal_or_morethan_5tpm_merged.txt"), sep="\t", header=True, index=True)


""" Perform 2 sample KS test to see if the distribution b/w 2 samples is significantly different:"""
tf_bound_tpm_array = df_prom_rnaseq_gene_merged_tf_bound["Value"].as_matrix()
tf_unbound_tpm_array = df_prom_rnaseq_gene_merged_tf_unbound["Value"].as_matrix()

stats.ks_2samp(tf_bound_tpm_array, tf_unbound_tpm_array)
ks_pval = stats.ks_2samp(tf_bound_tpm_array, tf_unbound_tpm_array)[1]


""" Concatenate 2 dataframes with TF bound and TF unbound to get them on single platforms """
df_prom_rnaseq_gene_merged_tf_bound["anno"] = "TF bound"
df_prom_rnaseq_gene_merged_tf_unbound["anno"] = "TF unbound"

df_tf_bound_unbound_combined = pd.concat([df_prom_rnaseq_gene_merged_tf_bound, df_prom_rnaseq_gene_merged_tf_unbound], ignore_index=True)
df_tf_bound_unbound_combined.to_csv(join(output_dir, "All_TF_bound_unbound_promoter_nearby_genes_w_tpm_vals.txt"), sep="\t", header=True, index=True)
