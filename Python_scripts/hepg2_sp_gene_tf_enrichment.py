import pandas as pd, numpy as np
import pybedtools, pickle, scipy 
import os,re,glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from os.path import join
from os.path import basename
from os.path import splitext

# Download RNAseq data from protein atlas:
output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/liver_hepg2_gene_analysis"
rna_tissue= "/Users/suryachhetri/Dropbox/for_genemodels/liver_hepg2_gene_analysis/rna_tissue.tsv"
rna_celline="/Users/suryachhetri/Dropbox/for_genemodels/liver_hepg2_gene_analysis/rna_celline.tsv"

tissue_df = pd.read_csv(rna_tissue, sep="\t"); tissue_df.head()
celline_df = pd.read_csv(rna_celline, sep="\t");celline_df.head()

liver_tissue_df = tissue_df[tissue_df["Sample"] == "liver"];liver_tissue_df.head()
hepg2_celline_df = celline_df[celline_df["Sample"] == "Hep G2"];hepg2_celline_df.head()

liver_hepg2_df = pd.merge(hepg2_celline_df, liver_tissue_df, on=["Gene", "Gene name"])
liver_hepg2_df_1 = liver_hepg2_df[(liver_hepg2_df["Value_x"] >= 1) & (liver_hepg2_df["Value_y"]>=1)]
liver_hepg2_df_1.rename(columns={"Gene name":"Gene_name"}, inplace=True)
liver_hepg2_df_1 = liver_hepg2_df_1.reset_index(drop=True)
liver_hepg2_df_1.to_excel(join(output_dir, "liver_hepg2_df.xls"), header=True, sheet_name="Sheet1")


null_sequence_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/unique_hepg2_analysis_total"
# main_dir = os.path.expanduser("~/Dropbox/for_chris/batch_I/hepg2_liver_gene_sp_analysis")

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/hepg2_liver_gene_sp_analysis"
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

output_dir =  join(main_dir, "files_2kb_from_all_tss_gt_1TPM")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

refgene_file = "/gpfs/gpfs1/home/schhetri/for_chris/refGene_hg19"
# refgene_file = os.path.expanduser("~/Dropbox/for_chris/refGene_hg19")

""" All TF containing dir """
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
# dir_path = os.path.expanduser("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/idr_passed_peaks_total/unique_TFs/SL*")
all_tf_file_list = glob.glob(dir_path) # represents all tf file list


""" Read liver and hepg2 specific genes """
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/hepg2_liver_gene_sp_analysis/liver_hepg2_df.xls"
# excelfile = os.path.expanduser("~/Dropbox/for_chris/Liver_specific_genes_HepG2.xlsx")
df_xls = pd.ExcelFile(excelfile).parse()
sp_gene_list = df_xls["Gene_name"].tolist()


def final_refgene_coords(refgene_hg19_file, sp_gene_list=False):
    refgene_hg19_file = refgene_file
    refgene_df = pd.read_csv(refgene_hg19_file, sep="\t")
    select_cols = ["name", "name2", "chrom", "strand", "txStart", "txEnd"]
    select_df = refgene_df.loc[:,select_cols]
    chrom_list = ["chr" + str(num) for num in range(1,23)] + ["chrX"] + ["chrY"]
    select_df_final = select_df[select_df["chrom"].isin(chrom_list)]

    """ use hepg2 specific gene list to extract the coordinates from refgene hg19 """
    if sp_gene_list:
        final_gene_df = select_df_final[select_df_final["name2"].isin(sp_gene_list)]
    else:
        final_gene_df = select_df_final

    final_gene_list = final_gene_df["name2"].unique().tolist()
    print "Total unique genes in your current refgene model is:", len(final_gene_list)

    final_norm_df =  final_gene_df.groupby(["chrom","strand","name2"]).apply(lambda x: (x["txStart"].min(), x["txEnd"].max()))
    final_norm_df =  final_norm_df.reset_index()
    final_norm_df[["chrom_start", "chrom_end"]] = final_norm_df.loc[:,0].apply(pd.Series)
    final_norm_df = final_norm_df.drop([0], axis=1)
    final_norm_df.sort_values(["name2"])
    final_refgene_coords = final_norm_df.drop_duplicates()

    return(final_refgene_coords)

refgene_coord_df = final_refgene_coords(refgene_file) 
refgene_coord_df = final_refgene_coords(refgene_file, sp_gene_list) 
final_sp_gene_list = refgene_coord_df["name2"].unique()


def generate_tss_updownstream_coords(refgene_coordinates_info, upstream_range, downstream_range):
    #upstream = 1000
    #downstream = 1000
    #upstream_range = upstream
    #downstream_range = downstream

    refgene_coordinates_info = refgene_coord_df
    final_refgene_df =  refgene_coordinates_info.sort_values(["chrom","chrom_start","chrom_end"])
    
    """ For positive strand, i.e if strand1 == "+": 
    tss_midpoint = final_refgene_df["chrom_start1"] """

    refgene_pos_df = final_refgene_df[final_refgene_df["strand"] == "+"]
    refgene_pos_df["tss_midpoint"] = final_refgene_df["chrom_start"]
    
    refgene_pos_df["chrom_start"] = refgene_pos_df["tss_midpoint"] + (-upstream_range)
    refgene_pos_df["chrom_end"] = refgene_pos_df["tss_midpoint"] + (downstream_range)

    """ For negative strand, i.e if strand1 == "-": 
    tss_midpoint = final_refgene_df["chrom_end"]
    Meth_r concept, start and end switched, so as to always maintain the higher coords for end site """

    refgene_neg_df = final_refgene_df[final_refgene_df["strand"] == "-"]
    refgene_neg_df["tss_midpoint"] = final_refgene_df["chrom_end"]
    
    # refgene_neg_df["chrom_start"] = refgene_neg_df["tss_midpoint"] - (-upstream_range) # originally correct
    # refgene_neg_df["chrom_end"] = refgene_neg_df["tss_midpoint"] - (downstream_range) # originally correct 
    refgene_neg_df["chrom_start"] = refgene_neg_df["tss_midpoint"] - (downstream_range) #but to maintain start as lower coord
    refgene_neg_df["chrom_end"] = refgene_neg_df["tss_midpoint"] - (-upstream_range) #but to maintain end as higher coord

    ### Combine the positive and negative stranded genes:
    refgene_model_df = pd.concat([refgene_pos_df, refgene_neg_df])
    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'name2', u'tss_midpoint', u'strand']
    tss_coord_df = refgene_model_df.loc[:,select_cols]
    tss_coord_df.to_csv(join(output_dir, "liver_hepg2_sp_gene_tss_coordinate_info_2kb.bed"), sep="\t", index=False, header=False)

    return(tss_coord_df)

tss_coord_df = generate_tss_updownstream_coords(refgene_coord_df, 2000, 2000)


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
        each_uniq_df = tss_coord_df[tss_coord_df["name2"].isin([each_gene])]
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

with open(join(output_dir,"final_master_gene_dict_all_genes_1TPM.pkl"), "w") as outfile:
    pickle.dump(master_gene_dict, outfile)

### To load back the pickled python objects:
# with open(join(output_dir,"final_master_gene_dictas_dict_3tf.pkl")) as infile:
#   df = pickle.load(infile)
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
final_ideas_df.to_csv(join(output_dir,"final_tf_ideas_peak_count_with_bg.bed"), sep="\t", header=True, index=False)

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
#final_pval_qval_df["percent_overlap"] = (final_pval_qval_df["sig_hits"]/final_pval_qval_df["total_peaks"])*100
#sig_hits_sorted_df = final_pval_qval_df.sort_values("sig_hits", ascending=False)

# Log transformation:
final_pval_qval_df["-log2(pvalue)"] = -np.log2(final_pval_qval_df["pvalue"].astype(float))
final_pval_qval_df["-log2(pval_corr)"] = -np.log2(final_pval_qval_df["pval_corr"].astype(float))
final_pval_qval_df_sorted = final_pval_qval_df.replace(-0, 0)

with open(join(output_dir,"final_pval_qval_df.pkl"), "w") as outfile:
    pickle.dump(final_pval_qval_df_sorted, outfile)

with open(join(output_dir,"final_pval_qval_df.pkl")) as infile:
    final_pval_qval_df_sorted = pickle.load(infile)

""" Enrichmnent value heatmap """
final_pval_qval_df_sorted.to_csv(join(output_dir,"final_tf_gene_enrichment_with_pval_qval.bed"), sep="\t", header=True, index=False)
final_pval_qval_heatmap_01 = final_pval_qval_df_sorted.pivot(index="tf_name", columns="gene_name", values = "-log2(pvalue)")
final_pval_qval_heatmap_01.to_csv(join(output_dir,"final_tf_gene_enrichment_pval_corr_heatmap_data.bed"), sep="\t", header=True, index=True)

# Based on the output of heatmap decided to exclude some of the datas:
sig_gene_list = final_pval_qval_df[final_pval_qval_df["pvalue"] < 0.05]["gene_name"].unique().tolist()
final_pval_qval_gene_filtered = final_pval_qval_df_sorted[final_pval_qval_df_sorted["gene_name"].isin(sig_gene_list)]
final_pval_qval_gene_filtered.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_gene_filtered.bed"), sep="\t", header=True, index=True)
final_pval_qval_heatmap = final_pval_qval_gene_filtered.pivot(index="tf_name", columns="gene_name", values = "-log2(pvalue)")
final_pval_qval_heatmap.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_gene_filtered.bed"), sep="\t", header=True, index=True)


################################################

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_factor_count_with_tss_gt_1TPM"
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

