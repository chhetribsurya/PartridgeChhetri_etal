import os, errno
import re
from glob import glob
from os.path import splitext, join, basename, dirname, expanduser
import pandas as pd, numpy as np
import pybedtools


# tf_files = expanduser("~/Dropbox/for_genemodels/idr_passed_peaks_total/test_analysis/SL*narrowPeak*")
tf_files = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"

# genecode_file = expanduser("~/Dropbox/for_genemodels/gencode.v19.annotation.gtf")
genecode_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/for_genecodeTssProfile/gencode.v19.annotation.gtf"

# processing genecode(GTF) file:
genecode_file_edited = join(dirname(genecode_file), "gencode.v19.annotation_edit1.gtf")

if not os.path.exists(genecode_file_edited):
	os.system( ''' tail -n+6 %s | awk 'BEGIN{print "chrom\tstart\tend\tgene_name\tgene_id\tstrand\tfeature"} \
			{if($3=="gene")print $1,$4,$5,$18,$10,$7,$3}' OFS="\t" > %s ''' %(genecode_file, genecode_file_edited))

# output_dir = os.path.expanduser("~/Dropbox/for_genemodels")
output_dir = os.path.expanduser("~/for_chris/batch_I/for_genecodeTssProfile/for_genemodels")
suboutput_dir = join(output_dir, "tss_tf_intersect_files")

# download HepG2 FPKM file from encode on fly:
download_dir = output_dir
if not os.path.exists(join(download_dir,"ENCFF139ZPW.tsv")):
	os.system("wget -P {} https://www.encodeproject.org/files/ENCFF139ZPW/@@download/ENCFF139ZPW.tsv".format(download_dir))

# Useful for job submission - if multiple threads running with a race condition to create the dir:
try:
	os.makedirs(suboutput_dir)
except OSError as exc:
	if exc.errno != errno.EEXIST:
		raise
	pass


def final_genecode_model(genecode_edited_file):
	""" 
	Function to filter raw genecode model for strand specific coordinate settings 
	
	Args:
		data_raw(TSV bedfile): Raw genecodeV19 gtf formatted bed file with selective cols
	
	Returns:
		data(pd.DataFame): Processed and cleansed genecode model data
	"""

	# Read Tab-separated value bedfiles:
	df = pd.read_csv(genecode_file_edited, sep="\t")
	regex_pat = re.compile(r";", flags=0) # flags=re.IGNORECASE
	df["gene_name"] = df["gene_name"].str.replace(regex_pat, "")
	df["gene_id"] = df["gene_id"].str.replace(regex_pat, "")
	df_pos = df.loc[df["strand"] == "+"]; df_pos.head()
	df_neg = df.loc[df["strand"] == "-"]; df_neg.head()
	df_neg_new = df_neg.loc[ :,["chrom","end","start", "gene_name", \
							"gene_id", "strand", "feature"]]; df_neg_new.head()
	df_neg_new.columns = df_pos.columns
	df_final = pd.concat([df_pos, df_neg_new]).sort_values(["chrom", "start","end"])
	print("Current dimension of the genecode model:{}".format(df_final.shape))
	final_genecode_model = df_final.drop_duplicates()
	print("Dropping duplicates,if any - current dimension of genecode model:{}\n\n".format(df_final.shape))
	final_genecode_model.to_csv(join(output_dir, "GenecodeV19_gtf_model_final.bed"), sep="\t", header = True, index= False)

	# If interested in tss only, just retain TSS (start site) with stop (=start+1); 
	# as pybed gets upset with (start > stop) for neg stranded genes
	final_tss_onlygenecode_model = final_genecode_model.copy()
	final_tss_onlygenecode_model["end"] = final_tss_onlygenecode_model["start"]+1
	final_tss_onlygenecode_model["geneID"] = final_tss_onlygenecode_model["gene_id"].str.replace(r"\..*", ""); final_tss_onlygenecode_model.head()
	final_genecode_model.to_csv(join(output_dir, "GenecodeV19_gtf_model_final.bed"), sep="\t", header = True, index= False)
	return(final_genecode_model, final_tss_onlygenecode_model)


genecode_coordinate_df, genecode_tssOnly_df = final_genecode_model(genecode_file_edited)


def final_gene_expression_model(gene_expression_file):
	""" 
	Function to filter raw gene expression file for downstream analysis 
	Args:
		data_raw(TSV bedfile): Raw gene_expression ENCODE bed file with selective cols
	Returns:
		data(pd.DataFame): Processed and cleansed gene expression model data
	"""

	expression_df = pd.read_csv(gene_expression_file, sep="\t")
	expression_df["geneID_ensemble"] = expression_df["gene_id"].str.replace(r"\..*", ""); expression_df.head()
	
	return(expression_df)

gene_exp_df = final_gene_expression_model(join(download_dir,"ENCFF139ZPW.tsv")); gene_exp_df.head(700).tail(5)


def find_nearest_tssDist_nearest_geneExp(feature_element_or_peakfile, gencode_tss_coords_df, gene_expression_df, custom=True):
	#gencode_tss_coords_df = genecode_tssOnly_df
	#gene_expression_df = gene_exp_df
	# Note: Set custom = True (if general bed file is used instead of Hotsite bed); else set False
	# Read peak bedfile
	peak_df = pd.read_csv(feature_element_or_peakfile, sep="\t")
	new_header = peak_df.columns.tolist()[:3]
	old_header = peak_df.columns.tolist()[3:]
	new_header = ["chrom", "chromStart", "chromEnd"]
	updated_header = new_header + old_header
	peak_df.columns = updated_header
	
	# For midpoint based peak file
	peak_df["mid"] = ((peak_df["chromStart"] +  peak_df["chromEnd"])/2).astype(int)
	peak_df["chromStart"] = peak_df["mid"]
	peak_df["chromEnd"] = peak_df["mid"] + 1
	peak_df.drop(["mid"], axis=1, inplace=True)
	
	peak_df_sorted = peak_df.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
	peak_header = peak_df_sorted.columns
	peak_pybed = pybedtools.BedTool.from_dataframe(peak_df_sorted)

	# Formatting of gencode tss df and find closest tss dist to feature element:
	gencode_tss_sorted_df = gencode_tss_coords_df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
	gencode_tss_header = gencode_tss_sorted_df.add_suffix('_gencode').columns
	gencode_tss_pybed = pybedtools.BedTool.from_dataframe(gencode_tss_sorted_df)
	peak_genecode_tss_closest = peak_pybed.closest(gencode_tss_pybed, D="b", t="first")
	peak_closest_tss_df = pd.read_csv(peak_genecode_tss_closest.fn, sep="\t", header=None)
	
	final_header = peak_header.tolist() + gencode_tss_header.tolist() + ["tss_dist"]
	peak_closest_tss_df.columns = final_header
	select_cols = ["chrom", "chromStart", "chromEnd", "uniq_tfcount", "binned_tf_count_1", "gene_name_gencode", "geneID_gencode", "tss_dist"]
	peak_tssfinal_df = peak_closest_tss_df.loc[:,select_cols]

	# Merge features carrying closest TssDist and gene with RNAseq gene expression file:
	merged_peak_tssdist_genexp = pd.merge(peak_tssfinal_df, gene_expression_df[["geneID_ensemble", "TPM", "FPKM"]], left_on = "geneID_gencode", right_on="geneID_ensemble", how="left")

	if custom:
		peak_tssfinal_df = peak_closest_tss_df.iloc[:,[0,1,2,5,-6,-2,-1]]
		peak_tssfinal_df.columns = ["chrom", "chromStart", "chromEnd", "uniq_tfcount", "gene_name_gencode", "geneID_gencode", "tss_dist"]
		# Merge features carrying closest TssDist and gene with RNAseq gene expression file:
		merged_peak_tssdist_genexp = pd.merge(peak_tssfinal_df, gene_expression_df[["geneID_ensemble", "TPM", "FPKM"]], left_on = "geneID_gencode", right_on="geneID_ensemble", how="left")
		return(merged_peak_tssdist_genexp)

	return(merged_peak_tssdist_genexp)

feature_element_or_peakfile = join(output_dir, "MergedBedswCounts.txt")
peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(feature_element_or_peakfile, genecode_tssOnly_df, gene_exp_df, custom=True)
peak_nearest_tss_genexp.to_csv(join(output_dir, "Hot_NonHotSites_nearest_tss_gene_expr_rynes.txt"), sep="\t", header=True, index=False)

# feature_element_or_peakfile = "/home/surya/Dropbox/for_genemodels/motifs_compinfo_analysis/MergedBedswCounts.txt"
feature_element_or_peakfile = join(output_dir, "Hotmotif_sites.bed")
peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(feature_element_or_peakfile, genecode_tssOnly_df, gene_exp_df, custom=False)
peak_nearest_tss_genexp.to_csv(join(output_dir, "Hot_NonHotSites_nearest_tss_gene_expr.txt"), sep="\t", header=True, index=False)


from plotnine import *
import pandas as pd
from os.path import join

""" Local machine plotting """

plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hot_NonHotSites_nearest_tss_gene_expr_rynes.txt", sep="\t")
plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hot_NonHotSites_nearest_tss_gene_expr.txt", sep="\t")

# Binning if needed :
bins_1 = [1,5,10,20,30,40,50,70,100,500]
names_1 = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
plot_df['binned_tf_count_1'] = pd.cut(plot_df["uniq_tfcount"], bins_1, right=False, labels=names_1)

# Order the factor/categorical variable to color legend accordingly:
plot_df["cobound_tf_bins_new"] = pd.Categorical(plot_df["binned_tf_count_1"], categories=["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"], ordered=True)
plot_df["log10FPKM"] = np.log10(plot_df["FPKM"]+1)
out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(y="log10FPKM", x="cobound_tf_bins_new", fill="cobound_tf_bins_new")+ 
        geom_boxplot(stat='boxplot',  outlier_shape="None") +
        ggtitle("Nearest gene expression") +
        theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Gene expression log10(FPKM)") + xlab("Number of TFs co-bound") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="TFs Co-Bound")) + 
        theme(
        axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot_title = element_text(size=14, face="bold"),
        legend_title = element_text(size=8, face="bold") 
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        ) 
    )
plot

ggsave(plot,join(out_dir, "Genexpresssion_cobound_sites_motifs.pdf"))

plot_df["cobound_tf_bins_new"] = pd.Categorical(plot_df["binned_tf_count_1"], categories=["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"], ordered=True)
plot_df["log10tss_dist"] = np.log10(abs(plot_df["tss_dist"]+1))
out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(y="log10tss_dist", x="cobound_tf_bins_new", fill="cobound_tf_bins_new")+ 
        geom_boxplot(stat='boxplot',  outlier_shape="None") +
        ggtitle("Nearest TSS dist") +
        theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("TSS distance log10(TssDist)") + xlab("Number of TFs co-bound") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="TFs Co-Bound")) +
        theme(
        axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot_title = element_text(size=14, face="bold"),
        legend_title = element_text(size=8, face="bold")
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
    )

plot

ggsave(plot,join(out_dir, "TssDist_cobound_sites_motifs.pdf"))
