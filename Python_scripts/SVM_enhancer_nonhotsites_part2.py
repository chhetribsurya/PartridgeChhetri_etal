import re, os 
from os.path import join, basename, expanduser
from glob import glob
from scipy import stats
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists


# IDEAS segmentation file:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")

# Peak based Hotsites (Rynes):
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm_peaklevel"
output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm_peaklevel"
final_hotpeak_df = pd.read_csv(join(output_dir, "MergedBedswCounts.txt"), sep="\t")
final_hotpeak_df.columns = ["chrom", "chromStart", "chromEnd", "total_peaks", "TFs_bound", "num_tfbound"]
final_hotpeak_df_70 = final_hotpeak_df.loc[final_hotpeak_df["num_tfbound"] > 70]
select_cols = ["chrom", "chromStart", "chromEnd", "num_tfbound"]
final_hotpeak_df_70 = final_hotpeak_df_70.loc[:, select_cols]
final_hotpeak_df_70.to_csv(join(output_dir, "Master_MergedBedswCounts_70more.txt" ), sep="\t", header=None, index=False)


def assign_IDEAS_State(peaks_motif_df):
    ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
    ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")
    peaks_motif_df = peaks_motif_df.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)

    ### Intersect Hepg2 bed file and Filter the rows with max base overlap size for duplicated rows:
    hepg2_ideas_df = pd.read_csv(ideas_hepg2_file, sep="\t")
    hepg2_ideas_df = hepg2_ideas_df.iloc[:, 1:5]
    hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "chromStart"])

    hepg2_bed_file = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
    peaks_motif_bed_file = pybedtools.BedTool.from_dataframe(peaks_motif_df)
    motif_ideas_intersect = peaks_motif_bed_file.intersect(hepg2_bed_file, wao=True)

    ### Note or remove - if needed -  the feature with 0 overlap size(last column)
    motif_ideas_intersect_df = pd.read_csv(motif_ideas_intersect.fn, sep="\t", header=None)
    intersect_df = motif_ideas_intersect_df.copy()
    # intersect_df =  motif_ideas_intersect_df[motif_ideas_intersect_df.iloc[:, -1] > 0] 
    
    ### Filter the rows with max base overlap size for duplicated rows:(column 8 here represents the size for base overlap)
    duplicated_df = intersect_df[intersect_df.duplicated([0,1,2,3], keep=False)]
    last_col_ID = duplicated_df.columns.tolist()[-1]
    duplicated_filtered_df = duplicated_df.loc[duplicated_df.groupby([0,1,2,3])[last_col_ID].idxmax()]
    non_duplicated_df = intersect_df[~intersect_df.duplicated([0,1,2,3], keep=False)] # non-duplicated df
    if (intersect_df.shape[0] == duplicated_df.shape[0] + non_duplicated_df.shape[0]): # 47743 = 34578 + 13165 
        print "Duplicated Rows filtered after groupby and merged successfully"

    motif_uniq_df = pd.concat([duplicated_filtered_df, non_duplicated_df], ignore_index=True)
    second_last_col_ID = duplicated_df.columns.tolist()[-2]
    final_col_select = range(0,len(peaks_motif_df.columns)) + [second_last_col_ID] #range(0,4)
    motif_ideas_final_df = motif_uniq_df.iloc[:,final_col_select]
    motif_ideas_final_df.columns = peaks_motif_df.columns.tolist() + ["ideas_state"]
    motif_ideas_final_df = motif_ideas_final_df.sort_values(["chrom", "chromStart", "chromEnd"])
    motif_ideas_final_df = motif_ideas_final_df.reset_index(drop=True) # reordering of the index
    motif_ideas_distribution = motif_ideas_final_df["ideas_state"].value_counts()

    # motif_ideas_final_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.bed"), header=True, index=True, sep="\t")
    # motif_ideas_final_df.to_pickle(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.pkl"))
    # motif_ideas_distribution.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_motifsfraction_ideas_piechart_dist.txt"), header=True, index=True, sep="\t")
    return(motif_ideas_final_df)


# Assign IDEAS annotation to sorted and all tf combined fimo coords:
final_hotpeak_df = final_hotpeak_df.iloc[:,[0,1,2,5]]
final_hotpeak_df = assign_IDEAS_State(final_hotpeak_df)

### For whole genome annotation:
replace_ideas = {"Tss" : "Prom assoc", "TssF" :  "Prom assoc", 'TssCtcf': "Prom assoc", 'PromCtcf': "Prom assoc", 
"TssW" : "Prom assoc", "PromP": "Prom assoc", 'PromF1': "Prom assoc", 'PromF2': "Prom assoc",
"Pol2" : "Gene body", "Gen5" : "Gene body", "Gen3": "Gene body", "Gen3Ctcf" : "Gene body", "Elon" : "Gene body", "ElonW" : "Gene body",
"CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", 
"Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh",
"Repr1" : "Heterochrom_repressed", "Repr2" : "Heterochrom_repressed", "ReprD":  "Heterochrom_repressed", "LowReprW":  "Heterochrom_repressed", 
"Low1":  "Heterochrom_repressed", "Low2":  "Heterochrom_repressed", "Quies" :  "Heterochrom_repressed", "Art1":  "Heterochrom_repressed",
"Art2":  "Heterochrom_repressed", "Zero" : "Heterochrom_repressed", 
"DnaseD1": "Euchromatin assoc", "DnaseD2": "Euchromatin assoc", "FaireW1": "Euchromatin assoc", "FaireW2": "Euchromatin assoc"}
remap_ideas = {"Prom assoc" : "Promoter associated", "Weak Enh" : "Enhancer associated", "Strong Enh" : "Enhancer associated"}

final_hotpeak_df["ideas_anno"] = final_hotpeak_df["ideas_state"].replace(replace_ideas)
final_hotpeak_df["ideas_reanno"] = final_hotpeak_df["ideas_anno"].replace(remap_ideas)
final_hotpeak_df.to_csv(join(output_dir, "final_hotpeak_df_with_Ideas.txt"), sep="\t", index=False, header=True)
final_hotpeak_df.to_pickle(join(output_dir, "final_hotpeak_df_with_Ideas.pkl"))

# Direct read from pickle:
final_hotpeak_df = pd.read_pickle(join(output_dir,"final_hotpeak_df_with_Ideas.pkl"))
final_hotpeak_df_enh = final_hotpeak_df[final_hotpeak_df["ideas_anno"] == "Strong Enh"]

# For centrimo enrichment run and scoring hotsites:
# Random 5676 enhancer non-hotsites (final_hotpeak_df.loc[(final_hotpeak_df["num_tfbound"]>70)].shape)
final_hotpeak_df_enh_samp = final_hotpeak_df_enh.loc[(final_hotpeak_df_enh["num_tfbound"]>=2)&(final_hotpeak_df_enh["num_tfbound"]<=10)]
final_hotpeak_df_enh_samp = final_hotpeak_df_enh_samp.sample(5676, random_state=20)
final_hotpeak_df_enh_samp.to_csv(join(output_dir, "final_hotpeak_df_enhancer_nonhotsites_2_10TF.txt"), sep="\t", index=False, header=False)
final_hotpeak_df_enh_samp_bed = final_hotpeak_df_enh_samp.iloc[:,0:3]
final_hotpeak_df_enh_samp_bed = final_hotpeak_df_enh_samp_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
final_hotpeak_df_enh_samp_bed.to_csv(join(output_dir, "NonhotsitesEnh_2_10TF.bed"), sep="\t", index=False, header=False)

# Random 5676 enhancer whole genome:
ideas_hepg2_df = pd.read_csv("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome", sep="\t")
ideas_hepg2_df.rename(columns={"name":"ideas_state"}, inplace=True)
ideas_hepg2_df["ideas_anno"] = ideas_hepg2_df["ideas_state"].replace(replace_ideas)
ideas_hepg2_df["ideas_reanno"] = ideas_hepg2_df["ideas_anno"].replace(remap_ideas)
ideas_hepg2_df = ideas_hepg2_df.iloc[:,1:]
ideas_hepg2_df.to_csv(join(output_dir, "ideas_hepg2_whole_genome.txt"), sep="\t", index=False, header=True)
ideas_hepg2_df.to_pickle(join(output_dir, "ideas_hepg2_whole_genome.pkl"))

ideas_hepg2_df_enh = ideas_hepg2_df[ideas_hepg2_df["ideas_anno"] == "Strong Enh"]
ideas_hepg2_df_enh_samp = ideas_hepg2_df_enh.sample(5676, random_state=20)
ideas_hepg2_df_enh_samp.to_csv(join(output_dir, "Hepg2_enhancer_whole_genome.txt"), sep="\t", index=False, header=False)
ideas_hepg2_df_enh_samp_bed = ideas_hepg2_df_enh_samp.iloc[:,0:3]
ideas_hepg2_df_enh_samp_bed = ideas_hepg2_df_enh_samp_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
ideas_hepg2_df_enh_samp_bed.to_csv(join(output_dir, "RandomGenomeEnh.bed"), sep="\t", index=False, header=False)

# Hotsites with 70 or more tf bound (n=5676) 
# final_hotpeak_df_hot = final_hotpeak_df_hot.sample(5676, random_state=20)
final_hotpeak_df_hot = final_hotpeak_df.loc[(final_hotpeak_df["num_tfbound"]>70)]
final_hotpeak_df_hot.to_csv(join(output_dir, "final_hotpeak_df_enhancer_hotsites_70_moreTF.txt"), sep="\t", index=False, header=False)
final_hotpeak_df_hot_bed = final_hotpeak_df_hot.iloc[:,0:3]
final_hotpeak_df_hot_bed = final_hotpeak_df_hot_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
final_hotpeak_df_hot_bed.to_csv(join(output_dir, "HotsitesEnh_70_moreTF.bed"), sep="\t", index=False, header=False)

# Mean SVM score for each 208 factor:
SVM_prediction_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm_peaklevel/hot_nonhot_random_scores"
SVM_prediction_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm_peaklevel/hot_nonhot_random_scores"
hot_score_filelist = glob(join(SVM_prediction_dir, "*.hot.scores"))
nonhot_score_filelist = glob(join(SVM_prediction_dir, "*.nonhot.scores"))
random_score_filelist = glob(join(SVM_prediction_dir, "*.random.scores"))
hot_score_filedict = {re.split(r"_gkmpredict.*.scores", basename(each))[0] : each for each in hot_score_filelist}
nonhot_score_filedict = {re.split(r"_gkmpredict.*.scores", basename(each))[0] : each for each in nonhot_score_filelist}
random_score_filedict = {re.split(r"_gkmpredict.*.scores", basename(each))[0] : each for each in random_score_filelist}

# SVM scores for each tf (self based score):
tf_based_svm_score_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm_peaklevel/random5000_samples/gkm_predict_output"
tf_based_svm_score_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm_peaklevel/gkm_predict_output"
tf_based_svm_filelist = glob(join(tf_based_svm_score_dir, "*.gkmpredict.scores"))
tf_based_score_filedict = {re.split(r".gkmpredict.scores", basename(each))[0] : each for each in tf_based_svm_filelist}


def compute_mean_svm_score(score_file_list, anno):
    tf_name_list = []
    mean_svm_score = []
    # score_file_list = hot_score_filelist
    for each in score_file_list:
        tf_name = re.split(r"_gkmpredict.*.scores", basename(each))[0]
        df = pd.read_csv(each, sep="\t", header=None)
        df.columns = ["loci", "score"]
        mean_score = df["score"].mean()
        mean_svm_score.append(mean_score)
        tf_name_list.append(tf_name)

    svm_score_df = pd.DataFrame({"tf_name": tf_name_list, anno + "_mean_score" : mean_svm_score})
    select_cols = ["tf_name", anno + "_mean_score"]
    svm_score_df = svm_score_df.loc[:,select_cols]
    return(svm_score_df)

hot_score_df = compute_mean_svm_score(hot_score_filelist, "hot")
nonhot_score_df = compute_mean_svm_score(nonhot_score_filelist, "nonhot")
random_score_df = compute_mean_svm_score(random_score_filelist, "random")
df_list = [hot_score_df, nonhot_score_df, random_score_df]

merged_df = reduce(lambda left, right : pd.merge(left, right, on=["tf_name"]), df_list )
merged_df.to_csv(join(output_dir, "hot_nonhot_random_scores_boxplot_data.txt"), sep="\t", header=True, index=False)


## Local machine plotting:

library(data.table)
library(dplyr)
library(ggplot2)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment"
read_df <- fread(file.path(output_dir, "hot_nonhot_random_scores_boxplot_data.txt"), sep="\t")
data_matrix <- read_df[,2:ncol(read_df)] %>% as.matrix()

ks_test_hot_nonhot = ks.test(read_df$hot_mean_score, read_df$nonhot_mean_score) # p-value = 5.966e-11
ks_test_hot_random = ks.test(read_df$hot_mean_score, read_df$random_mean_score) # p-value < 2.2e-16

pdf(file.path(output_dir, "hot_nonhot_random_scores_boxplot.pdf"))
boxplot(data_matrix, xaxt="n", col=c("red", "green", "blue"), main="Mean SVM Score Distribution", ylab="Mean SVM Score for sites(n=5676)")
legend("topright", legend = c("HOTsites", "NonHOTsites", "Random Enh"), col=c("red", "green", "blue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,3,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)", "Random Genome") , tick=FALSE , cex=0.3)
# legend("bottom", legend = c("HOTsites", "NonHOTsites", "Random Enh"), inset=c(0,1), xpd=TRUE, col=c("red", "green", "blue"), pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = T)
dev.off()

# Precision-recall base plot:
plot_df <- fread("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/PRAUC_all_TF_annotated_peaklevel.txt", sep="\t")
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment"
pdf(file.path(output_dir, "Mean_PRAUC_peaklevel_10X_nullsets.pdf"))
boxplot(plot_df$pr_auc ~ plot_df$Category, main="Mean PR-AUC=0.66", col=c("blue", "green"))
dev.off()


def merge_factors_svm_score(score_file_list):
    tf_svm_score_df = []
    for each in score_file_list:
        tf_name = re.split(r"_gkmpredict.*.scores", basename(each))[0]
        df = pd.read_csv(each, sep="\t", header=None)
        df.columns = ["loci", tf_name]
        tf_svm_score_df.append(df)

    merged_score_df = reduce(lambda left, right : pd.merge(left, right, on=["loci"]), tf_svm_score_df )
    return(merged_score_df)

combined_hot_score_df = merge_factors_svm_score(hot_score_filelist)
combined_nonhot_score_df = merge_factors_svm_score(nonhot_score_filelist)
combined_random_score_df = merge_factors_svm_score(random_score_filelist)


def compute_svm_score_rank(factors_score_dataframe):
    # For hotsite problem 2 # Rank hotsites with TFs_svmscore or classifier value or less : ()
    factors_score_dataframe = combined_hot_score_df
    df_rank = factors_score_dataframe.iloc[:,1:].rank(ascending=False, method="dense", pct=True)

    # Total frequency of hotsites containing top 5% classifier value for any TF present:
    df_rank["5perc_present"] = df_rank.apply(lambda row : row <= 0.05, axis=0).astype(int).apply(np.sum, axis=1)
    df_rank["75perc_present"] = df_rank.apply(lambda row : row > 0.25, axis=0).astype(int).apply(np.sum, axis=1)

    df_top5 = df_rank.reset_index(drop=True)[["5perc_present"]]
    df_top5["percent_classifier"] = "top5_percent"
    df_top5.rename(columns={"5perc_present": "num_bound_tfs"}, inplace=True)

    df_bottom75 = df_rank[["75perc_present"]]
    df_bottom75["percent_classifier"] = "bottom75_percent"
    df_bottom75.rename(columns={"75perc_present": "num_bound_tfs"}, inplace=True)
    df_rank_top_bottom = pd.concat([df_top5, df_bottom75], ignore_index=True)
    return(df_rank_top_bottom)

combined_hot_score_rank_df = compute_svm_score_rank(combined_hot_score_df)
combined_nonhot_score_rank_df = compute_svm_score_rank(combined_nonhot_score_df)
combined_random_score_rank_df = compute_svm_score_rank(combined_random_score_df)

combined_hot_score_rank_df.to_csv(join(output_dir,"Hotmotif_sites_histogram_data.txt"), header=True,  index=False, sep="\t")
combined_nonhot_score_rank_df.to_csv(join(output_dir,"NonHotmotif_sites_histogram_data.txt"), header=True,  index=False, sep="\t")
combined_random_score_rank_df.to_csv(join(output_dir,"Random_sites_histogram_data.txt"), header=True,  index=False, sep="\t")


def compute_tf_hot_svm_score_rank(hot_score_filedict, tf_based_score_filedict):
    tf_svm_score_df = []
    for tf, filepath in tf_based_score_filedict.iteritems():
        tf_name = re.split(r".gkmpredict.scores", basename(filepath))[0]
        tf_df = pd.read_csv(filepath, sep="\t", header=None)
        tf_df.columns = ["loci", tf_name]
        
        # Accessing hot_based tf svm scores:
        hot_tf_file = hot_score_filedict[tf_name]
        hot_df = pd.read_csv(hot_tf_file, sep="\t", header=None)
        hot_df.columns = ["loci", tf_name]
        hot_df["rank"] = map(lambda x: stats.percentileofscore(tf_df[tf_name], x), hot_df[tf_name])
        hot_df[tf_name + ".rank"] = abs(hot_df["rank"] - 100)/100 # since larger SVM should have lower rank percentile
        select_cols = ["loci", tf_name + ".rank"]
        hot_df = hot_df.loc[:,select_cols]
        tf_svm_score_df.append(hot_df)
        print("Processed scores for {}\n".format(tf_name))

    #map(lambda x: stats.percentileofscore([1, 2, 3, 4,5,6,7], x), [3,3,3])
    #map(lambda x: stats.percentileofscore(df["tf_scores"], x), df["tf_hotsite_svm_score"])
    df_rank = reduce(lambda left, right : pd.merge(left, right, on=["loci"]), tf_svm_score_df )

    df_rank["10perc_present"] = df_rank.iloc[:,1:].apply(lambda row : row <= 0.10, axis=0).astype(int).apply(np.sum, axis=1)
    df_rank["75perc_present"] = df_rank.iloc[:,1:].apply(lambda row : row > 0.25, axis=0).astype(int).apply(np.sum, axis=1)

    df_top5 = df_rank.reset_index(drop=True)[["10perc_present"]]
    df_top5["percent_classifier"] = "top10_percent"
    df_top5.rename(columns={"10perc_present": "num_bound_tfs"}, inplace=True)

    df_bottom75 = df_rank[["75perc_present"]]
    df_bottom75["percent_classifier"] = "bottom75_percent"
    df_bottom75.rename(columns={"75perc_present": "num_bound_tfs"}, inplace=True)
    df_rank_top_bottom = pd.concat([df_top5, df_bottom75], ignore_index=True)

    return(df_rank_top_bottom)

combined_hot_score_rank_df_1 = compute_tf_hot_svm_score_rank(hot_score_filedict, tf_based_score_filedict)
combined_hot_score_rank_df_1.to_csv(join(output_dir,"Hotmotif_sites_histogram_data_tfbased_10.txt"), header=True,  index=False, sep="\t")

combined_nonhot_score_rank_df_1 = compute_tf_hot_svm_score_rank(nonhot_score_filedict, tf_based_score_filedict)
combined_nonhot_score_rank_df_1.to_csv(join(output_dir,"NonHotmotif_sites_histogram_data_tfbased_10.txt"), header=True,  index=False, sep="\t")

combined_random_score_rank_df_1 = compute_tf_hot_svm_score_rank(random_score_filedict, tf_based_score_filedict)
combined_random_score_rank_df_1.to_csv(join(output_dir,"Random_sites_histogram_data_tfbased_10.txt"), header=True,  index=False, sep="\t")

###########################################################################
###########################################################################

# Find identifier(tf_name)from CisBP for all HepG2 factors expressed 1 TPM or more:
infile = open("/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme", "rb")
# infile = open("/gpfs/gpfs1/home/schhetri/Tools/final_custom_CisBP_database/Homosapiens_custom_cisbp_allmotifs.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        print(line)

cisbp_df = pd.Series(motif_list).str.split(expand=True)
cisbp_df.columns = ["motif", "cisbp_motif_id", "tf_name"]
cisbp_df["tf_name"] = cisbp_df["tf_name"].str.replace(r"_.*", "").str.replace(r"\(|\)", "").str.upper()

# All expressed HepG2 genes with 1 TPM or more:
rna_celline="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment/rna_celline.tsv"
celline_df = pd.read_csv(rna_celline, sep="\t");celline_df.head()
hepg2_celline_df = celline_df[celline_df["Sample"] == "Hep G2"];hepg2_celline_df.head()
hepg2_1tpm_df = hepg2_celline_df[hepg2_celline_df["Value"] >= 1]
hepg2_1tpm_df.rename(columns={"Gene name":"Gene_name"}, inplace=True)
hepg2_1tpm_df["Sample"] = hepg2_1tpm_df["Sample"].str.replace(" ", "")
hepg2_1tpm_df.to_csv(join(output_dir, "HepG2_factors_and_genes_with_1tpm_more.txt"), sep="\t", header=False, index=False)


# Annotate cisbp motifid from gene name (protein atlas) by merge with cisbp df:
hepg2_tfexpr_cisbp_df = pd.merge(hepg2_1tpm_df,cisbp_df, left_on="Gene_name", right_on="tf_name")
hepg2_tfexpr_cisbp_df_srt = hepg2_tfexpr_cisbp_df.sort_values(["Value", "tf_name"], ascending=False).reset_index(drop=True)
hepg2_tfexpr_cisbp_df_srt.to_csv(join(sub_outdir, "HepG2_exprGenes_Cisbp_factors_merged.txt"), sep="\t", header=True, index=False)
hepg2_tfexpr_cisbp_df_srt["cisbp_motif_id"].to_csv(join(sub_outdir, "HepG2_exprGenes_Cisbp_factors_with_Motif_Identifiers.txt"), header=False, index=False)


#############################

# For kmer stuff:
final_hotmotif_df_enh_samp = final_hotmotif_df_enh.loc[(final_hotmotif_df_enh["num_tfbound"]>=6)&(final_hotmotif_df_enh["num_tfbound"]<=10)]
final_hotmotif_df_enh_samp = final_hotmotif_df_enh_samp.sample(5000, random_state=20)
final_hotmotif_df_enh_samp.to_csv(join(sub_outdir, "final_hotmotif_df_enhancer_nonhotistes_2_10TF.txt"), sep="\t", index=False, header=False)
final_hotmotif_df_enh_samp = final_hotmotif_df_enh_samp.loc[:, ["num_tfbound", "hotsite_idx"]]

# filename = join(output_dir, "kmer_svm", "cv_nullseq_svmscores.pkl")    
# with open(filename, "rb") as readobj:
#     cv_svmscoredict = pickle.load(readobj)

filename = join(output_dir, "kmer_svm", "tf_svmscores.pkl")    
with open(filename, "rb") as readobj:
   master_tf_svmscoredict = pickle.load(readobj)

# For hotsite problem 1:
num_tfbound_list = []
tfscore_list = []
tf_namelist = []


# For hotsite problem 2:
master_list = []

# Hotsites dataframe - for problem 1 and 2 both:
for idx,row in final_hotmotif_df_enh_samp.iterrows():
    tfid_splitted=row["hotsite_idx"].split(",")
    num_tfbound =row["num_tfbound"]
    for each_id in tfid_splitted:
        # print num_tfbound, each_id
        tf_name = each_id.split("|")[0]
        if tf_name in master_tf_svmscoredict:
            if master_tf_svmscoredict[tf_name].get(each_id):
                tf_svmscore = master_tf_svmscoredict[tf_name].get(each_id)["score"]
                tfscore_list.append(tf_svmscore)
                num_tfbound_list.append(num_tfbound)
                tf_namelist.append(tf_name)
                # For hotsite problem 2:
                master_list.append([idx, num_tfbound, tf_svmscore, each_id, tf_name])
            else:
                tfscore_list.append(None)
                num_tfbound_list.append(num_tfbound)
                tf_namelist.append(tf_name)
                # For hotsite problem 2:
                master_list.append([idx, num_tfbound, None, each_id, tf_name])

# For hotsite problem 1:
tf_svmscore_df = pd.concat([pd.Series(num_tfbound_list), pd.Series(tfscore_list), pd.Series(tf_namelist)], axis=1)
tf_svmscore_df.columns = ["tf_cobound", "svm_score", "tf_name"]

# Binning HOT motif-sites
bins = [1,2,3,4,5,10,20,30,40,50,70,100,500]
names = ["1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
bins_1 = [1,5,10,20,30,40,50,70,100,500]
names_1 = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
tf_svmscore_df['binned_tf_count'] = pd.cut(tf_svmscore_df["tf_cobound"], bins, right=False, labels=names)
tf_svmscore_df['binned_tf_count_1'] = pd.cut(tf_svmscore_df["tf_cobound"], bins_1, right=False, labels=names_1)
# tf_svmscore_df.to_csv(join(output_dir,"Hotmotif_sites_problem1_data.txt"), header=True,  index=False, sep="\t")


