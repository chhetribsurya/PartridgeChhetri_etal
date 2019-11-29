
#######################################
# HOT motif Sites SVM train:         ##
#######################################

## Generate 5000 random sample for each fimo-bed file 

import re, os, pickle
from os import makedirs
from os.path import exists, join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"

suboutput_dir = "kmer_svm/unique_TFs_motif"
if not os.path.exists(join(output_dir, suboutput_dir)):
    os.makedirs(join(output_dir, suboutput_dir))

random_train_dir = "kmer_svm/random5000_samples"
if not os.path.exists(join(output_dir, random_train_dir)):
    os.makedirs(join(output_dir, random_train_dir))

# Generate dict for fimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_filelist = glob(file_pat)

# Generate dict for peak file list:
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"
bed_filelist = glob(file_pat)
bed_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_filedict = {bed_regex_pat.findall(file)[0]:file for file in bed_filelist}


###### Random sampling at motif level:
def svm_fimo_motifs_model(fimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)
    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["tf_name"] = tf_name
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)

    # Increases search-space up-and-downstream of motifs:
    # df["chromStart"] = df["chromStart"] -2
    # df["chromEnd"] = df["chromEnd"] +2
    df["midpos"] = ((df["chromStart"] + df["chromEnd"])/2).astype(int)
    df["chromStart"] = df["midpos"] - 15
    df["chromEnd"] = df["midpos"] + 15
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_id", "tf_name", "strand"]
    motif_select_df = df.loc[:,select_cols]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)
    return(motif_sorted_df)

regex_pat = re.compile(r'.*fimo_motif_(.*).txt$')

#fimo_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL1060_SE_VS_SL1167_fimo_motif_TCF12.txt"
for fimo_file in fimo_filelist: #
    TF_name = regex_pat.findall(basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    bed_file = bed_filedict[TF_name]
    print("Currently processing : {} TF\n".format(TF_name))
    fimo_coord_df = svm_fimo_motifs_model(fimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)     
    fimo_coord_df.to_csv(join(output_dir, suboutput_dir, TF_name + "_motifs.bed"), sep="\t", header=False, index=False) 
    
    # Random sampling of 5000 motif-sites for SVM train:
    if len(fimo_coord_df) > 5000:
        np.random.seed(10)
        sample_range = np.arange(0, len(fimo_coord_df))
        rindex = np.random.choice(sample_range, 5000, replace=False) # random permutation
        fimo_randn_df = fimo_coord_df.loc[rindex]
        # Make sure header = False; as nullgenerate seq would be upset with header:
        fimo_randn_df.to_csv(join(output_dir, random_train_dir, TF_name + "_motifs_sample.bed"), sep="\t", header=False, index=False)   
    
    else:
        fimo_coord_df.to_csv(join(output_dir, random_train_dir, TF_name + "_motifs_sample.bed"), sep="\t", header=False, index=False)   


# Set output dirs:
suboutput_dir = "kmer_svm_peaklevel/unique_TFs_motif"
if not os.path.exists(join(output_dir, suboutput_dir)):
    os.makedirs(join(output_dir, suboutput_dir))

random_train_dir = "kmer_svm_peaklevel/random5000_samples"
if not os.path.exists(join(output_dir, random_train_dir)):
    os.makedirs(join(output_dir, random_train_dir))

###### Random sampling at peak level:
def svm_fullpeak_model(peak_file):
    tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(peak_file))[0]
    df_read = pd.read_csv(peak_file, header=None, sep="\t")
    df_read = df_read.iloc[:,[0,1,2]]
    df_read.columns = ["chr", "start", "end"]
    df_read["tf_name"] = tf_name
    select_cols = ["chr","start","end", "tf_name"]
    df_read = df_read.loc[:,select_cols]
    print df_read.shape
    df_read = df_read.sort_values(select_cols).reset_index(drop=True)
    return(df_read)

regex_pat = re.compile(r".*narrowPeak_(.*)$")
for peak_file in bed_filelist: #
    TF_name = regex_pat.findall(basename(peak_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    bed_file = bed_filedict[TF_name]
    print("Currently processing : {} TF\n".format(TF_name))
    peak_coord_df = svm_fullpeak_model(peak_file)     
    peak_coord_df.to_csv(join(output_dir, suboutput_dir, TF_name + "_peaks.bed"), sep="\t", header=False, index=False) 
    
    # Random sampling of 5000 motif-sites for SVM train:
    if len(peak_coord_df) > 5000:
        np.random.seed(10)
        sample_range = np.arange(0, len(peak_coord_df))
        rindex = np.random.choice(sample_range, 5000, replace=False) # random permutation
        fimo_randn_df = peak_coord_df.loc[rindex]
        # Make sure header = False; as nullgenerate seq would be upset with header:
        fimo_randn_df.to_csv(join(output_dir, random_train_dir, TF_name + "_peaks_sample.bed"), sep="\t", header=False, index=False)   
    
    else:
        peak_coord_df.to_csv(join(output_dir, random_train_dir, TF_name + "_peaks_sample.bed"), sep="\t", header=False, index=False)   


# Set output dirs:
suboutput_dir = "kmer_svm_centpeaklevel/unique_TFs_motif"
if not os.path.exists(join(output_dir, suboutput_dir)):
    os.makedirs(join(output_dir, suboutput_dir))

random_train_dir = "kmer_svm_centpeaklevel/random5000_samples"
if not os.path.exists(join(output_dir, random_train_dir)):
    os.makedirs(join(output_dir, random_train_dir))

###### Random sampling at peak level:
def svm_centpeak_model(peak_file):
    tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(peak_file))[0]
    df_read = pd.read_csv(peak_file, header=None, sep="\t")
    df_read = df_read.iloc[:,[0,1,2]]
    df_read.columns = ["chr", "start", "end"]
    df_read["midpos"] = ((df_read["start"] + df_read["end"])/2).astype(int)
    df_read["start"] = df_read["midpos"] - 50
    df_read["end"] = df_read["midpos"] + 50
    df_read["tf_name"] = tf_name
    select_cols = ["chr","start","end", "tf_name"]
    df_read = df_read.loc[:,select_cols]
    print df_read.shape
    df_read = df_read.sort_values(select_cols).reset_index(drop=True)
    return(df_read)


regex_pat = re.compile(r".*narrowPeak_(.*)$")
for peak_file in bed_filelist: #
    TF_name = regex_pat.findall(basename(peak_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    bed_file = bed_filedict[TF_name]
    print("Currently processing : {} TF\n".format(TF_name))
    peak_coord_df = svm_centpeak_model(peak_file)     
    peak_coord_df.to_csv(join(output_dir, suboutput_dir, TF_name + "_centpeaks.bed"), sep="\t", header=False, index=False) 
    
    # Random sampling of 5000 motif-sites for SVM train:
    if len(peak_coord_df) > 5000:
        np.random.seed(10)
        sample_range = np.arange(0, len(peak_coord_df))
        rindex = np.random.choice(sample_range, 5000, replace=False) # random permutation
        fimo_randn_df = peak_coord_df.loc[rindex]
        # Make sure header = False; as nullgenerate seq would be upset with header:
        fimo_randn_df.to_csv(join(output_dir, random_train_dir, TF_name + "_centpeaks_sample.bed"), sep="\t", header=False, index=False)   
    
    else:
        peak_coord_df.to_csv(join(output_dir, random_train_dir, TF_name + "_centpeaks_sample.bed"), sep="\t", header=False, index=False)   


#################################################################
## Train SVM on random5000 samples and find PR-AUC for each TFs:
#################################################################

import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
# from sklearn.metrics import average_precision_score
# from sklearn.metrics import roc_auc_score

tf_namelist = []
PR_AUClist = []

output_dir= "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/random5000_samples/gkm_train_output/*.cvpred.txt"
crossval_filelist = glob(file_pat)

regex_pat = re.compile(r'^(.*)_motifs_sample.cvpred.txt')
for each in crossval_filelist:
    tf_name = regex_pat.findall(basename(each))[0]
    cv_svmscore_df = pd.read_csv(each, sep="\t", header=None)
    cv_svmscore_df.columns = ["score_idx","score", "class", "crossfold_valset"]
    cv_svmscore_df["score"] = cv_svmscore_df["score"].astype(float).round(4)
    select_cols = ["score", "class"]
    cv_svmscore_df = cv_svmscore_df.loc[:,select_cols]

    # Assign labels and scores predicted by clf to compute PR_AUC:
    y_test = cv_svmscore_df["class"]
    y_scores = cv_svmscore_df["score"]
    precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
    pr_areacurve = auc(recall, precision) # pr_areacurve = average_precision_score(y_test, y_scores)
    tf_namelist.append(tf_name)
    PR_AUClist.append(round(pr_areacurve,4))
    print('Area Under PR Curve(AP): {0:0.4f}'.format(pr_areacurve))

    # Alt method: The Area under Precision-Recall curve 
    # Alt method: The Area under an ROC(Receiver Operateing Characteristic) curve
    # fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    # roc_areacurve = auc(fpr, tpr)
    # Shortcut method if output is continuous probabitly; else use yscores = RandomForestClassifier.predict_proba(Xtest)[:,1]
    # pr_areacurve = average_precision_score(y_test, y_scores) 
    # roc_areacurve = roc_auc_score(y_test, y_scores)

pr_auc_df = pd.DataFrame({"tf_name":tf_namelist, "pr_auc":PR_AUClist})
pr_auclist_mean = round(pr_auc_df["pr_auc"].mean(),2)
print('Mean Area Under PR Curve(AP) for TF_list: {}'.format(pr_auclist_mean))
pr_auc_df.to_csv(join(output_dir, "PRAUC_all_TF.txt"), sep="\t", header=True, index=False)

# Give annotation to all TFs:
anno_df = pd.read_csv(join(output_dir, "TFs_Annotation_file.txt"), sep="\t")
prauc_anno_df = pd.merge(pr_auc_df, anno_df, left_on="tf_name",right_on="Target", how="left")
prauc_anno_df.to_csv(join(output_dir, "PRAUC_all_TF_annotated.txt"), sep="\t", header=True, index=False)

from plotnine import *
import pandas as pd
from os.path import join
""" Local machine plotting """

plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/PRAUC_all_TF_annotated.txt", sep="\t")

# Boxplot:
out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(y="pr_auc", x="annotation", fill="annotation")+ 
        geom_boxplot(stat='boxplot') +
        ggtitle("PRAUC distribution") +
        # theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Precision Recall AUC") + xlab("TF Category") +
        scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        # guides(fill=guide_legend(title="IDEAS Anno"))
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

ggsave(plot,join(out_dir, "PRAUC_ideas_anno.pdf"))

plot = (ggplot(plot_df) + 
        aes(y="pr_auc", x="Category", fill="Category")+ 
        geom_boxplot(stat='boxplot') +
        ggtitle("Mean PR-AUC : 0.74") +
        # theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Precision Recall AUC") + xlab("TF Category") +
        scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        # guides(fill=guide_legend(title="IDEAS Anno"))
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

ggsave(plot,join(out_dir, "PRAUC_dbf_crf.pdf"))

 
##################################################################
## Generate HOT motif-sites for downstream and svm-score analysis:
## using svm_fimo_motifs_model() func above, for *motifs.bed files
##################################################################

import pandas as pd, numpy as np
import pybedtools, pickle
from glob import glob

# Generate dict for fimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_filelist = glob(file_pat)

suboutput_dir = "kmer_svm/fimo_motifs_for_hotanalysis"
if not os.path.exists(join(output_dir, suboutput_dir)):
    os.makedirs(join(output_dir, suboutput_dir))

def hotsite_fimo_motifs_model(fimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)
    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["tf_name"] = tf_name
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_id", "tf_name", "strand"]
    motif_select_df = df.loc[:,select_cols]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)
    return(motif_sorted_df)

regex_pat = re.compile(r'.*fimo_motif_(.*).txt$')
for fimo_file in fimo_filelist: #
    TF_name = regex_pat.findall(basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    print("Currently processing : {} TF\n".format(TF_name))
    fimo_coord_df = hotsite_fimo_motifs_model(fimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)     
    fimo_coord_df.to_csv(join(output_dir, suboutput_dir, TF_name + "_motifs.bed"), sep="\t", header=False, index=False) 

# Generate hot-sites with proccessed fimo motif files: 
file_pat = join(output_dir, suboutput_dir, "*motifs.bed")
fimo_motif_filelist = glob(file_pat)

prehot_filelist = []
for each in fimo_motif_filelist:
    prehot_df = pd.read_csv(each, sep="\t",header=None)
    prehot_df.columns = ["chrom", "chromStart", "chromEnd", "motif_id", "tf_name", "strand"]
    linenum_df = pd.Series(prehot_df.index.values).astype(str)
    prehot_df["id"] = prehot_df["tf_name"] + "|"  + linenum_df
    prehot_df["motif_tag"] = prehot_df["tf_name"] + "|" + prehot_df["motif_id"]
    prehot_df["motif_linetag"] = prehot_df["tf_name"] + "|" + prehot_df["motif_id"] + "|" + linenum_df
    prehot_filelist.append(prehot_df)

# combine prehot dataframe for hotsite generation:
combined_prehot_df = pd.concat(prehot_filelist, ignore_index=True)
sorted_prehot_df = combined_prehot_df.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)

prehot_pybed = pybedtools.BedTool.from_dataframe(sorted_prehot_df)
merge_hot_pybed = prehot_pybed.merge(c=[5,5,8,8,9,7], o=["count","count_distinct","count", "count_distinct", "collapse", "collapse"])
merge_hot_pybed_df = pd.read_csv(merge_hot_pybed.fn, sep="\t", header=None) 
merge_hot_pybed_df.columns = ["chrom", "chromStart", "chromEnd", "total_tfcount", "uniq_tfcount", "total_motif_count", "distinct_motif_count", "merged_hotmotif_id", "id"]

select_cols = ["chrom", "chromStart", "chromEnd", "total_tfcount", "uniq_tfcount", "distinct_motif_count", "merged_hotmotif_id", "id" ]
final_hotmotif_df = merge_hot_pybed_df.loc[:,select_cols]
final_hotmotif_df = final_hotmotif_df.sort_values(["uniq_tfcount"]).reset_index(drop=True)

# Binning HOT motif-sites
bins = [1,2,3,4,5,10,20,30,40,50,70,100,500]
names = ["1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
bins_1 = [1,5,10,20,30,40,50,70,100,500]
names_1 = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
final_hotmotif_df['binned_tf_count'] = pd.cut(final_hotmotif_df["uniq_tfcount"], bins, right=False, labels=names)
final_hotmotif_df['binned_tf_count_1'] = pd.cut(final_hotmotif_df["uniq_tfcount"], bins_1, right=False, labels=names_1)
final_hotmotif_df.to_csv(join(output_dir, "Hotmotif_sites.bed"), sep="\t", header=True, index=False)
final_hotmotif_df.to_pickle(join(output_dir,"Hotmotif_sites.pkl"))

# Frequency count:
hotsite_freq_count = final_hotmotif_df["uniq_tfcount"].value_counts().reset_index(name="site_count").rename(columns={"index":"uniq_tfcount"})
binned_hotsite_freq_count = final_hotmotif_df["binned_tf_count"].value_counts().reset_index(name="site_count").rename(columns={"index":"uniq_tfcount"})


from plotnine import *
import pandas as pd
from os.path import join

""" Local machine plotting """

test = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hotmotif_sites.bed", sep="\t")
test["diff"] = test.chromEnd - test.chromStart 
test["diff"].mean(), test["diff"].max(), test["diff"].min()
plot_df = test["binned_tf_count_1"].value_counts().reset_index(name="site_count").rename(columns={"index":"uniq_tfcount"})
plot_df["log2(site_count)"] = np.log2(plot_df["site_count"])
plot_df["site_count"].sum()

# Order the factor/categorical variable to color legend accordingly:
plot_df["uniq_tfcount_new"] = pd.Categorical(plot_df["uniq_tfcount"], categories=["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"], ordered=True)

out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(y="log2(site_count)", x="uniq_tfcount_new")+ 
        geom_bar(stat='identity') +
        ggtitle("Hotsites dist with number of TFs cobound") +
        theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Log2(Site Counts)") + xlab("Unique TFs co-bound") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        # guides(fill=guide_legend(title="TFs Co-Bound")) +
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

ggsave(plot,join(out_dir, "Hotmotif_sites_barplot_distribution_Figure.pdf"))


##############################################################################
## Generate dictionary for TF SVM-scores for each TFs for downstream analysis:
##############################################################################

import re, pickle
from glob import glob
from os.path import join, basename

# Create SVMscore dict for cross-fold validation sets(only on null seq scores)
cv_svmscoredict = {}
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/random5000_samples/gkm_train_output/*.cvpred.txt"
crossval_filelist = glob(file_pat)

regex_pat = re.compile(r'^(.*)_motifs_sample.cvpred.txt')
for each in crossval_filelist:
    tf_name = regex_pat.findall(basename(each))[0]
    cv_svmscore_df = pd.read_csv(each, sep="\t", header=None)
    cv_svmscore_df.columns = ["score_idx","score", "class", "crossfold_valset"]
    cv_svmscore_df["score"] = cv_svmscore_df["score"].astype(float).round(4)
    select_cols = ["score", "class"]
    cv_svmscore_df = cv_svmscore_df.loc[:,select_cols]

    # Select only null seq i.e class -1 for later usage
    cv_svmscore_df = cv_svmscore_df.loc[cv_svmscore_df["class"] == -1 ]
    svm_scoredict = {}
    tf_svm_scoredict = cv_svmscore_df.to_dict(orient="list")
    svm_scoredict[tf_name] = tf_svm_scoredict
    cv_svmscoredict.update(svm_scoredict)

filename = join(output_dir, "kmer_svm", "cv_nullseq_svmscores.pkl")    
fileobj = open(filename, 'wb')
pickle.dump(cv_svmscoredict, fileobj)
fileobj.close()

# with open(filename, "rb") as readobj:
#     cv_svmscoredict = pickle.load(readobj)

# Create SVMscore dict for scored DNA sequence:
master_tf_svmscoredict = {}
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/random5000_samples/gkm_predict_output/*gkmpredict.scores"
filelist = glob(file_pat)

regex_pat = re.compile(r'^(.*)_motifs.gkmpredict.scores')
for each in filelist:
    tf_name = regex_pat.findall(basename(each))[0]
    svmscore_df = pd.read_csv(each, sep="\t", header=None)
    svmscore_df.columns = ["score_idx","score"]
    svmscore_df["score"] = svmscore_df["score"].astype(float).round(4)
    svmscore_df["id"] = tf_name + "|" + svmscore_df.index.astype(str)
    select_cols = ["score", "id"]
    svmscore_df = svmscore_df.loc[:,select_cols]

    if svmscore_df.shape[0] > 6000:
        np.random.seed(10)
        sample_range = np.arange(0, len(svmscore_df))
        rindex = np.random.choice(sample_range, 5000, replace=False) # random permutation
        random_sample_set = set(rindex)
        total_df_set = set(svmscore_df.index)
        testset_leftidx = list(total_df_set - random_sample_set)
        test_sample_df = svmscore_df.loc[testset_leftidx] # testset_df.to_dict(orient="list")
        
        svm_scoredict = {}
        tf_svm_scoredict = test_sample_df.set_index(["id"]).to_dict("index")
        svm_scoredict[tf_name] = tf_svm_scoredict
        master_tf_svmscoredict.update(svm_scoredict)

    else:
        svm_scoredict = {}
        tf_svm_scoredict = svmscore_df.set_index(["id"]).to_dict("index")
        svm_scoredict[tf_name] = tf_svm_scoredict
        master_tf_svmscoredict.update(svm_scoredict)

filename = join(output_dir, "kmer_svm", "tf_svmscores.pkl")    
fileobj = open(filename, 'wb')
pickle.dump(master_tf_svmscoredict, fileobj)
fileobj.close()

# with open(filename, "rb") as readobj:
#    master_tf_svmscoredict = pickle.load(readobj)


#######################################################################################
## Downstream analysis using hotsites-tfbound_file cv-and-tf svmscore dictionary above:
#######################################################################################
import pandas as pd
import re, pickle
from glob import glob
from os.path import join, basename

output_dir= "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
final_hotmotif_df = pd.read_pickle(join(output_dir,"Hotmotif_sites.pkl")) 
final_hotmotif_df.rename(columns={"uniq_tfcount" : "num_tfbound", "id" : "hotsite_idx"}, inplace=True) 
final_hotmotif_df = final_hotmotif_df.loc[:, ["num_tfbound", "hotsite_idx"]]
# pd.read_csv(join(output_dir, "Hotmotif_sites.bed"), sep="\t")

filename = join(output_dir, "kmer_svm", "cv_nullseq_svmscores.pkl")    
with open(filename, "rb") as readobj:
    cv_svmscoredict = pickle.load(readobj)

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
for idx,row in final_hotmotif_df.iterrows():
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
tf_svmscore_df.to_csv(join(output_dir,"Hotmotif_sites_problem1_data.txt"), header=True,  index=False, sep="\t")

# Create boxplot dataframe for hotsites and distribution of clf values:
nullseq_scorelist = []
for each_tf in cv_svmscoredict:
    nullseq_svmscore_df = pd.DataFrame(cv_svmscoredict[each_tf])
    nullseq_scorelist.append(nullseq_svmscore_df)

combined_nullseq = pd.concat(nullseq_scorelist, ignore_index=True)
combined_nullseq["binned_tf_count_1"] = "Matched_null"
combined_nullseq_df = combined_nullseq.loc[:,["binned_tf_count_1", "score"]]
combined_nullseq_df.columns = ["cobound_tf_bins", "svm_score"]
tf_svmscore_boxdf = tf_svmscore_df.loc[:,["binned_tf_count_1", "svm_score"]]
tf_svmscore_boxdf.columns = ["cobound_tf_bins", "svm_score"]

boxplot_svmscore_df = pd.concat([tf_svmscore_boxdf,combined_nullseq_df], ignore_index=True)
boxplot_svmscore_df.to_csv(join(output_dir,"Hotmotif_sites_problem1_boxplot_data.txt"), header=True,  index=False, sep="\t")


################################################


# For hotsite problem 2 (master list df) and TF_bound >=50
master_tf_svmscore_df = pd.DataFrame(master_list)
master_tf_svmscore_df.columns = ["final_hotsite_idx", "tf_bound", "score", "id", "tf_name"]
master_tf_svmscore_df = master_tf_svmscore_df.loc[(master_tf_svmscore_df["tf_bound"]>=50)]

# Handles the case with hotsites containing more than 1 peaks or motifs 
# for same factor like FOXA3 and KDM2A at hotsite 4 and 6 in hepg2 hotsites:
df_grouped = master_tf_svmscore_df.groupby(["final_hotsite_idx", "tf_bound", "tf_name"])["score"].mean().reset_index(name="svmscore")
df_final = df_grouped.pivot(index="final_hotsite_idx", columns="tf_name", values="svmscore")

# For hotsite problem 2 # Rank hotsites with TFs_svmscore or classifier value or less : ()
df_rank = df_final.rank(ascending=False,method="dense",pct=True)

# Total frequency of hotsites containing top 5% classifier value for any TF present:
df_rank["5perc_present"] = df_rank.apply(lambda row : row <= 0.05, axis=0).astype(int).apply(np.sum, axis=1)
df_rank["75perc_present"] = df_rank.apply(lambda row : row > 0.25, axis=0).astype(int).apply(np.sum, axis=1)
# df_rank_merged = pd.merge(df_rank.reset_index(), master_tf_svmscore_df.loc[:,["final_hotsite_idx",  "tf_bound"]], on = "final_hotsite_idx")
# df_rank_final = df_rank_merged.drop_duplicates()

df_top5 = df_rank.reset_index(drop=True)[["5perc_present"]]
df_top5["percent_classifier"] = "top5_percent"
df_top5.rename(columns={"5perc_present": "num_bound_tfs"}, inplace=True)

df_bottom75 = df_rank[["75perc_present"]]
df_bottom75["percent_classifier"] = "bottom75_percent"
df_bottom75.rename(columns={"75perc_present": "num_bound_tfs"}, inplace=True)

df_rank_top_bottom = pd.concat([df_top5, df_bottom75], ignore_index=True)
df_rank_top_bottom.to_csv(join(output_dir,"Hotmotif_sites_problem2_histogram_data.txt"), header=True,  index=False, sep="\t")


###############################################


# For hotsite problem 3 (master list df) - for piechart:
df_rank_melt = pd.melt(df_rank.reset_index(), id_vars=['final_hotsite_idx'], var_name = ["tf_name"], value_vars=df_rank.columns[:-2].tolist())
max_classifier_val_idx = df_rank_melt.groupby(["final_hotsite_idx"])["value"].idxmin()
hotsite_piechart_df = df_rank_melt.loc[max_classifier_val_idx]
hotsite_piechart_final_df = hotsite_piechart_df["tf_name"].value_counts()
hotsite_piechart_final_df = hotsite_piechart_final_df.reset_index(name="hotsite_count")

# Though total hotsite is # hotsite_piechart_df.shape[0] i.e 2040, giving 1
hotsite_piechart_final_df["total_hotsite"] = hotsite_piechart_final_df[hotsite_piechart_final_df["hotsite_count"]>=1].shape[0] # 2040 
hotsite_piechart_final_df["hotsite_fraction_w_recurring_motif"] = ((hotsite_piechart_final_df["hotsite_count"])/2040)*100 #
hotsite_piechart_final_df["reshaped_percent"] = ((hotsite_piechart_final_df["hotsite_count"])/1000)*100 #
hotsite_piechart_final_df.to_csv(join(output_dir,"Hotmotif_sites_problem3_piechart_data.txt"), header=True,  index=False, sep="\t")


###################### if needed; else ignore this analysis  #########################

# Grouping by TF bound gives us the frequency of hotsites with TF classifier value 
# Give score of 1 if present that is (more than 1 TF with classifier value 0.05); 
df_rank["5perc_binary"] = np.where(df_rank["5perc_present"] > 0, 1, 0 )
df_rank_merged = pd.merge(df_rank.reset_index(), master_tf_svmscore_df.loc[:,["final_hotsite_idx",  "tf_bound"]], on = "final_hotsite_idx")
df_rank_final = df_rank_merged.drop_duplicates()
df_rank_final.groupby(["tf_bound"])["5perc_binary"].sum().reset_index(name="5perc_present")

###################### if needed; else ignore above analysis  #########################

from plotnine import *
import pandas as pd
from os.path import join

""" Local machine plotting """

plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hotmotif_sites_problem1_boxplot_data.txt", sep="\t")

# Order the factor/categorical variable to color legend accordingly:
plot_df["cobound_tf_bins_new"] = pd.Categorical(plot_df["cobound_tf_bins"], categories=["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+", "Matched_null"], ordered=True)

out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(y="svm_score", x="cobound_tf_bins_new", fill="cobound_tf_bins_new")+ 
        geom_boxplot(stat='boxplot',  outlier_shape="None") +
        ggtitle("SVM weights distribution") +
        theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("SVM classifier scores") + xlab("Number of TFs co-bound") +
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
# plot
ggsave(plot,join(out_dir, "Hotmotif_sites_problem1_boxplot_svm_clf_weights_Figure.pdf"))


plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hotmotif_sites_problem2_histogram_data.txt", sep="\t")
out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot = (ggplot(plot_df) + 
        aes(x="num_bound_tfs", fill="percent_classifier")+ 
        geom_histogram(stat ='bin', binwidth=1) +
        ggtitle("Ranked Classifier-Weights Distribution") +
        # theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Number of Hotsites(>=50 TFs cobound)") + xlab("Number of bound TFs with SVM classifier values (each site)") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="Ranked Classifier Value")) +
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

ggsave(plot,join(out_dir, "Hotmotif_sites_problem2_histogram_svm_clf_value_figure.pdf"))

