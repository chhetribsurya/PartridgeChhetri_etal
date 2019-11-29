import pandas as pd, numpy as np
import pybedtools
import pickle, scipy
import os, re
from glob import glob
from os.path import join, basename, splitext, dirname
from scipy.stats import ks_2samp

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/gatad2a_cobind_analysis"
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

""" All TF containing dir """
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob(dir_path) # represents all tf file list

gatad2a_peakpath = glob(dir_path + "GATAD2A*")[0]
gatad2a_peakpath_random = "/gpfs/gpfs1/home/schhetri/for_encode/spp_idr_peak_call_scripts/pseudo_double_rep/pe_spp_scripts_newcluster_k562/motif_pipeline_k562/random_peaks/GATAD2A_peaks_random"

foxa3_peakpath = glob(dir_path + "FOXA3*")[0]
foxa3_peakpath_random = "/gpfs/gpfs1/home/schhetri/for_encode/spp_idr_peak_call_scripts/pseudo_double_rep/pe_spp_scripts_newcluster_k562/motif_pipeline_k562/random_peaks/FOXA3_peaks_random"

gatad2a_cobind_list = ["ZNF219", "SMAD4", "RARA", "ARID5B", "FOXA3", "SOX13_iso1"]
gatad2a_cobind_peakspath = [ glob(dir_path + each + "*")[0] for each in gatad2a_cobind_list ]


def str_join(df, sep, *cols):
    from functools import reduce
    return reduce(lambda x, y: x.astype(str).str.cat(y.astype(str), sep=sep), [df[col] for col in cols])


def pybedfile_format(file_path, tf_name):
    peak_df = pd.read_csv(file_path, sep="\t", header=None)
    peak_df = peak_df.iloc[:, 0:3]
    peak_df.columns = ["chrom", "chromStart", "chromEnd"]
    peak_df["tf_name"] = tf_name
    peak_df = peak_df.sort_values(["chrom", "chromStart", "chromEnd"])
    bedfile = pybedtools.BedTool.from_dataframe(peak_df)
    return(bedfile)


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

def find_nearest_tssDist_nearest_geneExp(feature_element_or_peakfile, gencode_tss_coords_df, gene_expression_df, custom=True):
    # gencode_tss_coords_df = genecode_tssOnly_df
    # gene_expression_df = gene_exp_df
    # Note: Set custom = True (if general bed file is used instead of Hotsite bed); else set False
    # Read peak bedfile
    peak_df = pd.read_csv(feature_element_or_peakfile, sep="\t", header=None)
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
        peak_tssfinal_df = peak_closest_tss_df.iloc[:,[0,1,2,-6,-2,-1]]
        peak_tssfinal_df.columns = ["chrom", "chromStart", "chromEnd", "gene_name_gencode", "geneID_gencode", "tss_dist"]
        # Merge features carrying closest TssDist and gene with RNAseq gene expression file:
        merged_peak_tssdist_genexp = pd.merge(peak_tssfinal_df, gene_expression_df[["geneID_ensemble", "TPM", "FPKM"]], left_on = "geneID_gencode", right_on="geneID_ensemble", how="left")
        return(merged_peak_tssdist_genexp)
    return(merged_peak_tssdist_genexp)


# genecode_file = expanduser("~/Dropbox/for_genemodels/gencode.v19.annotation.gtf")
genecode_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/for_genecodeTssProfile/gencode.v19.annotation.gtf"

# processing genecode(GTF) file:
genecode_file_edited = join(dirname(genecode_file), "gencode.v19.annotation_edit1.gtf")

if not os.path.exists(genecode_file_edited):
    os.system( ''' tail -n+6 %s | awk 'BEGIN{print "chrom\tstart\tend\tgene_name\tgene_id\tstrand\tfeature"} \
            {if($3=="gene")print $1,$4,$5,$18,$10,$7,$3}' OFS="\t" > %s ''' %(genecode_file, genecode_file_edited))

# output_dir = os.path.expanduser("~/Dropbox/for_genemodels")
output_dir = os.path.expanduser("~/for_chris/batch_I/for_genecodeTssProfile/for_genemodels")
# suboutput_dir = join(output_dir, "tss_tf_intersect_files")

# download HepG2 FPKM file from encode on fly:
download_dir = output_dir
if not os.path.exists(join(download_dir,"ENCFF139ZPW.tsv")):
    os.system("wget -P {} https://www.encodeproject.org/files/ENCFF139ZPW/@@download/ENCFF139ZPW.tsv".format(download_dir))

gatad2a_bedfile = pybedfile_format(gatad2a_peakpath, "GATAD2A")
gatad2a_bedfile_random = pybedfile_format(gatad2a_peakpath_random, "GATAD2A_random")
foxa3_bedfile = pybedfile_format(foxa3_peakpath, "FOXA3")
foxa3_bedfile_random = pybedfile_format(foxa3_peakpath_random, "FOXA3_random")

foxa3_only = foxa3_bedfile.intersect(gatad2a_bedfile, v=True)
gatad2a_only = gatad2a_bedfile.intersect(foxa3_bedfile, v=True)
gatad2a_foxa3 = gatad2a_bedfile.intersect(foxa3_bedfile, u=True)

genecode_coordinate_df, genecode_tssOnly_df = final_genecode_model(genecode_file_edited)
gene_exp_df = final_gene_expression_model(join(download_dir,"ENCFF139ZPW.tsv")); gene_exp_df.head(700).tail(5)

foxa3_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(foxa3_bedfile.fn, genecode_tssOnly_df, gene_exp_df, custom=True)
gatad2a_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(gatad2a_bedfile.fn, genecode_tssOnly_df, gene_exp_df, custom=True)

foxa3_only_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(foxa3_only.fn, genecode_tssOnly_df, gene_exp_df, custom=True)
gatad2a_only_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(gatad2a_only.fn, genecode_tssOnly_df, gene_exp_df, custom=True)
gatad2a_foxa3_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(gatad2a_foxa3.fn, genecode_tssOnly_df, gene_exp_df, custom=True)
foxa3_random_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(foxa3_bedfile_random.fn, genecode_tssOnly_df, gene_exp_df, custom=True)
gatad2a_random_peak_nearest_tss_genexp = find_nearest_tssDist_nearest_geneExp(gatad2a_bedfile_random.fn, genecode_tssOnly_df, gene_exp_df, custom=True)


df_dict = {"foxa3_only" : foxa3_only_peak_nearest_tss_genexp, 
            #"gatad2a_only" :gatad2a_only_peak_nearest_tss_genexp, 
            "gatad2a_foxa3" : gatad2a_foxa3_peak_nearest_tss_genexp,
            "foxa3_random" : foxa3_random_peak_nearest_tss_genexp, 
            #"gatad2a_random" : gatad2a_random_peak_nearest_tss_genexp
            }

gatad2a_foxa3_peak_nearest_tss_genexp["TPM"].mean() #37.43318338098503
foxa3_only_peak_nearest_tss_genexp["TPM"].mean() #30.249388613630046
gatad2a_only_peak_nearest_tss_genexp["TPM"].mean() #30.05162444654126
foxa3_random_peak_nearest_tss_genexp["TPM"].mean() #9.892717790588724
gatad2a_random_peak_nearest_tss_genexp["TPM"].mean() #30.05162444654126

ks_2samp(gatad2a_foxa3_peak_nearest_tss_genexp["TPM"], foxa3_only_peak_nearest_tss_genexp["TPM"])
ks_2samp(gatad2a_foxa3_peak_nearest_tss_genexp["TPM"], foxa3_random_peak_nearest_tss_genexp["TPM"])
ks_2samp(gatad2a_foxa3_peak_nearest_tss_genexp["TPM"], gatad2a_only_peak_nearest_tss_genexp["TPM"])
ks_2samp(gatad2a_random_peak_nearest_tss_genexp["TPM"], foxa3_random_peak_nearest_tss_genexp["TPM"])


concat_list = []
for key,val in df_dict.iteritems():
    val["annotation"] = key
    print(val.head())
    concat_list.append(val)

combined_df = pd.concat(concat_list, ignore_index=True)
combined_df.to_csv(join(output_dir, "gatad2a_foxa3_nearest_tss_gene_expr.txt"), sep="\t", header=True, index=False)

foxa3_only_peak_nearest_tss_genexp.iloc[:, 0:3].to_csv(join(output_dir, "foxa3_only_peaks.txt"), sep="\t", header=False, index=False)

foxa3_random_peak_nearest_tss_genexp.iloc[:, 0:3].to_csv(join(output_dir, "foxa3_random_peaks.txt"), sep="\t", header=False, index=False)

gatad2a_foxa3_peak_nearest_tss_genexp.iloc[:, 0:3].to_csv(join(output_dir, "gatad2a_foxa3_peaks.txt"), sep="\t", header=False, index=False)

gatad2a_only_peak_nearest_tss_genexp.iloc[:, 0:3].to_csv(join(output_dir, "gatad2a_only_peaks.txt"), sep="\t", header=False, index=False)

gatad2a_random_peak_nearest_tss_genexp.iloc[:, 0:3].to_csv(join(output_dir, "gatad2a_random_peaks.txt"), sep="\t", header=False, index=False)

genecode_tssOnly_df.to_csv(join(output_dir, "genecode_genes.txt"), sep="\t", header=True, index=False)

gene_exp_df.to_csv(join(output_dir, "hepg2_genes_exp.txt"), sep="\t", header=True, index=False)


###############################################################

## In local machine: (using GREAT output for nearest TSS dist)

input_dir = "/Users/suryachhetri/Dropbox/for_genemodels/revision_2"
genecode_genes = pd.read_csv(join(input_dir, "genecode_genes.txt"), sep="\t")
hepg2_genes_exp = pd.read_csv(join(input_dir, "hepg2_genes_exp.txt"), sep="\t")
genecode_hepg2_gene_exp = pd.merge(genecode_genes[[u'chrom', u'start', u'end', u'gene_name', u'geneID']], 
                                    hepg2_genes_exp[["geneID_ensemble", "TPM", "FPKM"]], 
                                    left_on = "geneID", right_on="geneID_ensemble", how="left"
                                    )

gatad2a_foxa3_exp = pd.read_csv(join(input_dir, "gatad2a_foxa3_peaks_genes.txt"), skiprows=1, header=None, sep="\t")
gatad2a_foxa3_exp.columns = ["name", "gene_name"]
gatad2a_foxa3_exp["gene_name"] = gatad2a_foxa3_exp["gene_name"].str.replace(r"(/d+)", "", regex=True)
gatad2a_foxa3_exp["gene_name_list"] = gatad2a_foxa3_exp["gene_name"].str.split(" ")
gatad2a_foxa3_exp_df = gatad2a_foxa3_exp["gene_name_list"].apply(pd.Series).iloc[:,0:2]
gatad2a_foxa3_exp_df.columns = ["gene_name", "tss_dist"]

foxa3_only_exp = pd.read_csv(join(input_dir, "foxa3_only_peaks_genes.txt"), skiprows=1, header=None, sep="\t")
foxa3_only_exp.columns = ["name", "gene_name"]
foxa3_only_exp["gene_name"] = foxa3_only_exp["gene_name"].str.replace(r"(/d+)", "", regex=True)
foxa3_only_exp["gene_name_list"] = foxa3_only_exp["gene_name"].str.split(" ")
foxa3_only_exp_df = foxa3_only_exp["gene_name_list"].apply(pd.Series).iloc[:,0:2]
foxa3_only_exp_df.columns = ["gene_name", "tss_dist"]
# foxa3_only_exp_df = pd.DataFrame(foxa3_only_exp["gene_name_list"].values.tolist(), columns = ["gene_name", "tss_dist"])

gatad2a_only_exp = pd.read_csv(join(input_dir, "gatad2a_only_peaks_genes.txt"), skiprows=1, header=None, sep="\t")
gatad2a_only_exp.columns = ["name", "gene_name"]
gatad2a_only_exp["gene_name"] = gatad2a_only_exp["gene_name"].str.replace(r"(/d+)", "", regex=True)
gatad2a_only_exp["gene_name_list"] = gatad2a_only_exp["gene_name"].str.split(" ")
gatad2a_only_exp_df = gatad2a_only_exp["gene_name_list"].apply(pd.Series).iloc[:,0:2]
gatad2a_only_exp_df.columns = ["gene_name", "tss_dist"]

foxa3_random_exp = pd.read_csv(join(input_dir, "foxa3_random_peaks_genes.txt"), skiprows=1, header=None, sep="\t")
foxa3_random_exp.columns = ["name", "gene_name"]
foxa3_random_exp["gene_name"] = foxa3_random_exp["gene_name"].str.replace(r"(/d+)", "", regex=True)
foxa3_random_exp["gene_name_list"] = foxa3_random_exp["gene_name"].str.split(" ")
foxa3_random_exp_df = foxa3_random_exp["gene_name_list"].apply(pd.Series).iloc[:,0:2]
foxa3_random_exp_df.columns = ["gene_name", "tss_dist"]

gatad2a_random_exp = pd.read_csv(join(input_dir, "gatad2a_random_peaks_genes.txt"), skiprows=1, header=None, sep="\t")
gatad2a_random_exp.columns = ["name", "gene_name"]
gatad2a_random_exp["gene_name"] = gatad2a_random_exp["gene_name"].str.replace(r"(/d+)", "", regex=True)
gatad2a_random_exp["gene_name_list"] = gatad2a_random_exp["gene_name"].str.split(" ")
gatad2a_random_exp_df = gatad2a_random_exp["gene_name_list"].apply(pd.Series).iloc[:,0:2]
gatad2a_random_exp_df.columns = ["gene_name", "tss_dist"]

genecode_exp.head()
gatad2a_foxa3_exp.head()
foxa3_only_exp.head()
gatad2a_only_exp.head()
foxa3_random_exp.head()
gatad2a_random_exp.head()

gatad2a_foxa3_tpm = pd.merge(genecode_hepg2_gene_exp[["gene_name", "geneID", "TPM", "FPKM"]], gatad2a_foxa3_exp_df[["gene_name"]])
foxa3_only_tpm = pd.merge(genecode_hepg2_gene_exp[["gene_name", "geneID", "TPM", "FPKM"]], foxa3_only_exp_df[["gene_name"]])
gatad2a_only_tpm = pd.merge(genecode_hepg2_gene_exp[["gene_name", "geneID", "TPM", "FPKM"]], gatad2a_only_exp_df[["gene_name"]])
foxa3_random_tpm = pd.merge(genecode_hepg2_gene_exp[["gene_name", "geneID", "TPM", "FPKM"]], foxa3_random_exp_df[["gene_name"]])
gatad2a_random_tpm = pd.merge(genecode_hepg2_gene_exp[["gene_name", "geneID", "TPM", "FPKM"]], gatad2a_random_exp_df[["gene_name"]])

gatad2a_foxa3_tpm["TPM"].median()
foxa3_only_tpm["TPM"].median()
gatad2a_only_tpm["TPM"].median()
foxa3_random_tpm["TPM"].median()
gatad2a_random_tpm["TPM"].median()

ks_2samp(gatad2a_foxa3_tpm["TPM"], foxa3_only_tpm["TPM"])
ks_2samp(gatad2a_foxa3_tpm["TPM"], foxa3_random_tpm["TPM"])
ks_2samp(gatad2a_foxa3_tpm["TPM"], gatad2a_only_tpm["TPM"])
ks_2samp(gatad2a_foxa3_tpm["TPM"], gatad2a_only_tpm["TPM"])
ks_2samp(gatad2a_foxa3_tpm["TPM"], gatad2a_random_tpm["TPM"])

df_list = [ gatad2a_foxa3_tpm["TPM"], foxa3_only_tpm["TPM"], gatad2a_only_tpm["TPM"], 
    foxa3_random_tpm["TPM"], gatad2a_random_tpm["TPM"] ]

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/revision_2"
df_concat = pd.concat(df_list, axis=1)
df_concat.columns = ["gatad2a_foxa3", "foxa3_only", "gatad2a_only", "foxa3_random", "gatad2a_random"]
df_concat.to_csv(join(output_dir, "gatad2a_foxa3_random_combined_gene_exp.txt"), sep="\t", header=True, index=False)

