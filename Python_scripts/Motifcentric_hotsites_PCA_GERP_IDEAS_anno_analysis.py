import pandas as pd
import numpy as np
import os
from os.path import join, expanduser, basename
import pybedtools
import pickle

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")
output_dir = expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis")
hotmotif_file = join(output_dir, "Hotmotif_sites.bed")
loci_df = pd.read_csv(hotmotif_file, sep="\t"); loci_df.head()
select_cols = ["chrom", "chromStart", "chromEnd", "uniq_tfcount", "merged_hotmotif_id", "id"]
final_loci_df = loci_df.loc[:, select_cols]


def assign_IDEAS_State(peaks_motif_df, cols_to_retain):
    ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
    # ideas_hepg2_file = os.path.expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/HepG2_ideas_whole_genome")
    
    # file = join(output_dir, "Hotmotif_sites.bed")
    # peaks_motif_df = pd.read_csv(file, sep = "\t")    
    col_select = cols_to_retain
    peaks_motif_df = peaks_motif_df.iloc[:, col_select ]
    peaks_motif_df = peaks_motif_df.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)

    ### Intersect Hepg2 bed file and Filter the rows with max base overlap size for duplicated rows:
    hepg2_ideas_df = pd.read_csv(ideas_hepg2_file, sep="\t")
    hepg2_ideas_df = hepg2_ideas_df.iloc[:, 1:5]
    hepg2_ideas_df = hepg2_ideas_df.sort_values(["chrom", "chromStart", "chromEnd"])

    hepg2_bed_file = pybedtools.BedTool.from_dataframe(hepg2_ideas_df)
    peaks_motif_bed_file = pybedtools.BedTool.from_dataframe(peaks_motif_df)

    motif_ideas_intersect = peaks_motif_bed_file.intersect(hepg2_bed_file, wao=True)

    ### Note or remove - if needed -  the feature with 0 overlap size(last column)
    motif_ideas_intersect_df = pd.read_csv(motif_ideas_intersect.fn, sep="\t", header=None)
    intersect_df = motif_ideas_intersect_df.copy()
    #intersect_df =  motif_ideas_intersect_df[motif_ideas_intersect_df.iloc[:, -1] > 0] 
    
    ### Filter the rows with max base overlap size for duplicated rows:(column 8 here represents the size for base overlap)
    duplicated_df = intersect_df[intersect_df.duplicated([0,1,2,3,4,5], keep=False)]
    last_col_ID = duplicated_df.columns.tolist()[-1]
    duplicated_filtered_df = duplicated_df.loc[duplicated_df.groupby([0,1,2,3,4,5])[last_col_ID].idxmax()]
    non_duplicated_df = intersect_df[~intersect_df.duplicated([0,1,2,3,4,5], keep=False)] # non-duplicated df
    if (intersect_df.shape[0] == duplicated_df.shape[0] + non_duplicated_df.shape[0]): # 47743 = 34578 + 13165 
        print "Duplicated Rows filtered after groupby and merged successfully"

    motif_uniq_df = pd.concat([duplicated_filtered_df, non_duplicated_df], ignore_index=True)
    second_last_col_ID = duplicated_df.columns.tolist()[-2]
    final_col_select = col_select + [second_last_col_ID]
    motif_ideas_final_df = motif_uniq_df.iloc[:,final_col_select]
    motif_ideas_final_df.columns = peaks_motif_df.columns.tolist() + ["ideas_state"]
    motif_ideas_final_df = motif_ideas_final_df.sort_values(["chrom", "chromStart", "chromEnd"])
    motif_ideas_final_df = motif_ideas_final_df.reset_index(drop=True) # reordering of the index
    motif_ideas_distribution = motif_ideas_final_df["ideas_state"].value_counts()

    motif_ideas_final_df.to_csv(join(output_dir, "Hotmotif_sites_ideas.bed"), header=True, index=True, sep="\t")
    motif_ideas_final_df.to_pickle(join(output_dir, "Hotmotif_sites_ideas.pkl"))
    motif_ideas_distribution.to_csv(join(output_dir, "Hotmotif_sites_ideas_piechart_dist.txt"), header=True, index=True, sep="\t")

    return(motif_ideas_final_df)


# Assign IDEAS annotation to sorted and all tf combined fimo coords:
retain_cols = [0,1,2,3,4,5]
hotmotif_ideas_df = assign_IDEAS_State(final_loci_df, retain_cols)
hotmotif_ideas_df = hotmotif_ideas_df.sort_values(["uniq_tfcount"]).reset_index(drop=True)
# hotmotif_ideas_df1 = pd.read_pickle(join(output_dir, "Hotmotif_sites_ideas.pkl"))

# Filter hotsites with >=10:
hotsite_df = hotmotif_ideas_df.loc[hotmotif_ideas_df["uniq_tfcount"] >= 2].reset_index(drop=True)
hotsite_df["hotloci_id"] = hotsite_df[["chrom", "chromStart", "chromEnd"]].apply(lambda row: "_".join(row.astype(str)), axis=1)
select_cols = ["hotloci_id", "uniq_tfcount", "ideas_state", "merged_hotmotif_id", "id"]
new_hotsite_df = hotsite_df.loc[:, select_cols]


# only using merged_hotmotif_id so dropping other id:
include_id = "id"
drop_id = "merged_hotmotif_id"
final_hotsite_df = new_hotsite_df.drop([drop_id], axis=1)
final_hotsite_df[include_id] = final_hotsite_df[include_id].str.replace(r'\|\d+', "")
final_hotsite_df[include_id] = final_hotsite_df[include_id].str.replace(r'\|', "-")
final_hotsite_df[include_id] = final_hotsite_df[include_id].str.replace(",", "|")
pca_cols = ["hotloci_id", "uniq_tfcount", "ideas_state"]

# For binary call:
tf_dummies = final_hotsite_df[include_id].str.get_dummies().add_prefix("TF_")
pca_hotsite_df = pd.concat([final_hotsite_df.loc[:,pca_cols], tf_dummies], axis=1)

# For sum calls instead of binary 0/1:
# expand_str = final_hotsite_df[include_id].str.split('|', expand=True)
# tf_dummies = pd.get_dummies(expand_str.stack(dropna=False)).sum(level=0).add_prefix("TF_")
# pca_hotsite_df = pd.concat([final_hotsite_df.loc[:,pca_cols], tf_dummies], axis=1)

## pca_hotsite_df = final_hotsite_df.loc[:,pca_cols].join(tf_dummies)
## pd.get_dummies(final_hotsite_df.set_index(pca_cols)[include_id].str.split('|', expand=True)\
##                  .stack(dropna=False), prefix='TF').sum(level=0)

subset_df = pca_hotsite_df["hotloci_id"].str.split("_", expand=True)
subset_df.rename(columns={0:"chrom", 1: "chromStart", 2:"chromEnd"}, inplace=True)
final_pca_df = subset_df.join(pca_hotsite_df)

# Remapping of ideas states:
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
remap_ctcf = {"TssCtcf": "CTCF bound", "PromCtcf" : "CTCF bound", "CtcfO" : "CTCF bound", "Gen3Ctcf" : "CTCF bound", "Ctcf" : "CTCF bound"}

final_pca_df["ideas_state"] = final_pca_df["ideas_state"].replace(replace_ideas)
final_pca_df.to_pickle(join(output_dir, "final_pca_df_2.pkl"))
final_pca_df_10tf = final_pca_df.loc[final_pca_df["uniq_tfcount"] >=10]
final_pca_df.to_pickle(join(output_dir, "final_pca_df_10_binary.pkl"))
# final_pca_df["ideas_reanno"] = final_pca_df["ideas_anno"]

# for ctcf based analysis:
final_pca_df_ctcf = final_pca_df.copy()
final_pca_df_ctcf["ideas_ctcf"] = final_pca_df_ctcf["ideas_state"].map(remap_ctcf)
values = {'ideas_ctcf': "CTCF not bound"}
final_pca_df_ctcf.fillna(value=values, inplace=True)
final_pca_df_ctcf.dropna(subset=["ideas_ctcf"], inplace=True)
final_pca_df_ctcf["ideas_state"] = final_pca_df_ctcf["ideas_ctcf"]
final_pca_df_ctcf.drop(["ideas_ctcf"], axis=1, inplace=True)


#################################################################
# PCA of samples based on TF's bound :
#################################################################

import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler
import matplotlib; import matplotlib.pyplot as plt
import os; from os.path import join, expanduser
from ggplot import *
from plotnine import *

""" Run locally after pickle file transfer from remote """
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10_binary.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))

feature_df = final_pca_df.iloc[:, 6:]
target_df = final_pca_df["ideas_state"]

# Scaling/Standardizing features 
# x = scale(feature_df.values) # df.values(Separating out the features with numpy array conversion)
# x = StandardScaler().fit_transform(x)

# PCA analysis (PCA projection to 2D/3D):
pca = PCA() # PCA(n_components=3)
principal_comp = pca.fit_transform(feature_df.values) 
# principal_comp = pca.fit_transform(x) # X already normalized with fraction
print("original shape:   ", x.shape)
print("transformed shape:", principal_comp.shape) # eigen decomposed data(reduced to n dim)

# Transform to df:
principal_comp_df = pd.DataFrame(principal_comp)

# The amount of variance that each PC explains
comp = pca.components_ # accessing components/eigen vectors
var = pca.explained_variance_ratio_ # var = pca.explained_variance_ # accessing eigen values or magnitude of eigen vectors

pc1 = round(var[0]*100, 2)
pc2 = round(var[1]*100, 2)
pc3 = round(var[2]*100, 2)
pc4 = round(var[3]*100, 2)

print("Amount of variance expained by PC1, PC2, PC3, PC4 : {}, {}, {}, {} respectively\n".format(pc1,pc2,pc3,pc4))
print("Amount of variance expained by each PC :\n{}\n".format(var))
print("Sum of variance:{}".format(sum(var)))

cum_var = np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
print("Total variance expained (cumsum) :\n{}".format(cum_var))

# Including 3rd dataframe for Color reference:
final_df = pd.concat([ principal_comp_df, target_df], axis = 1)

# Rename for later convenience:
final_df.rename(columns={0:"Principal_component_1", 1:"Principal_component_2", 2: "Principal_component_3", 3: "Principal_component_4"}, inplace=True)

# Remap ideas annotation for final plot:
map_ideas = {"Prom assoc" : "Promoter", "Strong Enh" : "Strong Enhancer", "Weak Enh" : "Weak Enhancer", "Ctcf assoc" : "CTCF assoc"}
final_df["ideas_state"] = final_df["ideas_state"].map(map_ideas)
values = {'ideas_state': "Other"}
final_df.fillna(value=values, inplace=True)
final_df.dropna(subset=["ideas_state"], inplace=True)

# Order the factor/categorical variable to color legend accordingly:
final_df["ideas_state"] = pd.Categorical(final_df["ideas_state"], categories=["Other",  "CTCF assoc", "Promoter", "Weak Enhancer", "Strong Enhancer"], ordered=True)
select_cols = ["Principal_component_1", "Principal_component_2", "Principal_component_3", "Principal_component_4", "ideas_state"]
final_df = final_df.loc[:,select_cols]
map_color = {"CTCF assoc" : "purple", "Other" : "green", "Promoter" : "red", "Strong Enhancer" : "orange", "Weak Enhancer" : "yellow"}
final_df["color_label"] = final_df["ideas_state"].map(map_color)
final_df.to_csv(join(output_dir, "PCA_analysis_hotmotif_sites_2_wo_scale.txt"), sep="\t", header=True, index=False)

# Plot principal components:
plot = (ggplot(final_df) + 
        aes("Principal_component_1", "Principal_component_2", fill="ideas_state")+ 
        geom_point(size=1) +
        # geom_text(label = final_df["ideas_state"], size=6) +
        # geom_text(label=final_df["Label"], size=8)+
        ggtitle("PCA analysis") +
        # theme_bw() +
        scale_x_continuous(name="Principal Component 1 " + "(" + str(pc1) + "%)") + scale_y_continuous(name="Principal Component 2 " + "(" + str(pc2) + "%)") +
        # xlab("Principal Component 1") + ylab("Principal Component 2") +
        scale_fill_manual(name="IDEAS Anno", values=["black","blueviolet", "red", "yellow", "orange" ]) +
        guides(fill=guide_legend(title="IDEAS States")) +
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

ggsave(plot, join(output_dir, "PCA_analysis_hotmotif_sites_fig_test.pdf"), width=6, height=6)


# For CTCF based analysis:
# Plot principal components:
feature_df = final_pca_df_ctcf.iloc[:, 6:]
target_df = final_pca_df_ctcf["ideas_state"]

# Separating out the features and convert to numpy arrays:
x = feature_df.values 

# Scaling the values:
x = scale(x)

# PCA analysis (PCA projection to 2D/3D):
pca = PCA() # PCA(n_components=3)
principal_comp = pca.fit_transform(x) # X already normalized with fraction

# Transform to df:
principal_comp_df = pd.DataFrame(principal_comp)

# The amount of variance that each PC explains
comp = pca.components_ # accessing components/eigen vectors
var = pca.explained_variance_ratio_ # var = pca.explained_variance_ # accessing eigen values or magnitude of eigen vectors

pc1 = round(var[0]*100, 2)
pc2 = round(var[1]*100, 2)
pc3 = round(var[2]*100, 2)

print("Amount of variance expained by PC1, PC2, PC3 : {}, {}, {} respectively\n".format(pc1,pc2,pc3))
print("Amount of variance expained by each PC :\n{}\n".format(var))
print("Sum of variance:{}".format(sum(var)))

cum_var = np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
print("Total variance expained (cumsum) :\n{}".format(cum_var))

# Including 3rd dataframe for Color reference:
final_df_ctcf = pd.concat([ principal_comp_df, target_df], axis = 1)

# Rename for later convenience:
final_df_ctcf.rename(columns={0:"Principal_component_1", 1:"Principal_component_2", 2: "Principal_component_3"}, inplace=True)
final_df_ctcf.to_csv(join(output_dir, "PCA_analysis_hotmotif_sites_ctcf.txt"), sep="\t", header=True, index=False)

# Order the factor/categorical variable to color legend accordingly:
final_df_ctcf["ideas_state"] = pd.Categorical(final_df_ctcf["ideas_state"], categories=["CTCF bound", "CTCF not bound"], ordered=True)

plot = (ggplot(final_df_ctcf) + 
        aes("Principal_component_2", "Principal_component_3", fill="ideas_state")+ 
        geom_point(size=1) +
        # geom_text(label = final_df["ideas_state"], size=6) +
        # geom_text(label=final_df["Label"], size=8)+
        ggtitle("PCA analysis") +
        # theme_bw() +
        scale_x_continuous(name="Principal Component 2 " + "(" + str(pc2) + "%)") + scale_y_continuous(name="Principal Component 3 " + "(" + str(pc3) + "%)") +
        # xlab("Principal Component 1") + ylab("Principal Component 2") +
        scale_fill_manual(name="IDEAS Anno", values=["blueviolet", "green"]) +
        guides(fill=guide_legend(title="IDEAS States")) +
        theme(
        axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot_title = element_text(size=14, face="bold"),
        legend_title = element_text(size=8, face="bold")
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
    )

ggsave(plot, join(output_dir, "PCA_analysis_hotmotif_sites_pc2_pc3_ctcf_fig.pdf"), width=6, height=6)


###################################
# t-SNE (high dimensional visualization)
""" Running locally after pickle file transfer from remote """

###################################

# final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))

map_ideas = {"Prom assoc" : "Promoter", "Strong Enh" : "Strong Enhancer", "Weak Enh" : "Weak Enhancer", "Ctcf assoc" : "CTCF assoc"}
final_pca_df["ideas_state"] = final_pca_df["ideas_state"].map(map_ideas)
values = {'ideas_state': "Other"}
final_pca_df.fillna(value=values, inplace=True)
final_pca_df.dropna(subset=["ideas_state"], inplace=True)

# Randomly sample 
final_pca_df_sample = final_pca_df.groupby('ideas_state').apply(lambda x: x.sample(2000, random_state=20)).reset_index(drop=True)

# tsne plot:
feature_df = final_pca_df_sample.iloc[:, 6:]
target_df = final_pca_df_sample["ideas_state"]

# Separating out the features and convert to numpy arrays:
x = feature_df.values 

# Separating out the target or labels: 
y = target_df.values

# Scaling the values:
x = scale(x)

# PCA analysis (PCA projection to 2D/3D):
pca_50 = PCA(n_components=50) # PCA(n_components=3)
principal_comp_50 = pca_50.fit_transform(x) # X already normalized with fraction

print("original shape:   ", x.shape)
print("transformed shape:", principal_comp_50.shape) # eigen decomposed data(reduced to n dim)
# Transform to df:
principal_comp_df_50 = pd.DataFrame(principal_comp_50)
# The amount of variance that each PC explains
comp = pca_50.components_ # accessing components/eigen vectors
var = pca_50.explained_variance_ratio_ # var = pca.explained_variance_ # accessing eigen values or magnitude of eigen vectors
pc1 = round(var[0]*100, 2)
pc2 = round(var[1]*100, 2)
pc3 = round(var[2]*100, 2)
print("Amount of variance expained by PC1, PC2, PC3 : {}, {}, {} respectively\n".format(pc1,pc2,pc3))
print("Amount of variance expained by each PC :\n{}\n".format(var))
print("Sum of variance for 50 components :{}".format(sum(var)))
cum_var = np.cumsum(np.round(pca_50.explained_variance_ratio_, decimals=4)*100)
print("Total variance expained (cumsum) :\n{}".format(cum_var))


import time
from sklearn.manifold import TSNE
time_start = time.time()

# with pca reduction to save computation time and cpu memory exhaustion:
tsne = TSNE(n_components=3, verbose=1, perplexity=40, n_iter=500)
tsne_pca_results = tsne.fit_transform(principal_comp_50)
tsne_df = pd.DataFrame(tsne_pca_results)
# Including 3rd dataframe for Color reference:
final_df = pd.concat([tsne_df, target_df], axis = 1)
final_df.rename(columns={0:"tsne_pca_1", 1:"tsne_pca_2", 2:"tsne_pca_3"}, inplace=True)
final_df.to_csv(join(output_dir, "tsne_analysis_hotmotif_sites_2000_pca.txt"), sep="\t", header=True, index=False)
# Order the factor/categorical variable to color legend accordingly:
final_df["ideas_annp"] = pd.Categorical(final_df["ideas_anno"], categories=["Other", "CTCF assoc", "Weak Enhancer", "Strong Enhancer", "Promoter"], ordered=True)
print 't-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start)

# Plot tsne principal components:
plot = (ggplot(final_df) + 
        aes("tsne_pca_1", "tsne_pca_2", fill="ideas_reanno")+ 
        geom_point(size=1) +
        # geom_text(label = final_df["ideas_state"], size=6) +
        # geom_text(label=final_df["Label"], size=8)+
        ggtitle("tSNE dimensions colored by Hotmotif site PCA") +
        # theme_bw() +
        scale_x_continuous(name="x-tsne") + scale_y_continuous(name = "y-tsne" ) +
        # xlab("Principal Component 1") + ylab("Principal Component 2") +
        scale_fill_manual(name="IDEAS Anno", values=["black", "blueviolet", "yellow", "orange", "red"]) +
        guides(fill=guide_legend(title="IDEAS States")) +
        theme(
        axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot_title = element_text(size=14, face="bold"),
        legend_title = element_text(size=8, face="bold")
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
    )

ggsave(plot, join(output_dir, "tsne_analysis_hotmotif_sites_2000_pca_fig.pdf"), width=7, height=7)


#########################################

# with out pca reduction - using original plot:
hotsite_df = final_pca_df[final_pca_df["uniq_tfcount"] >=50]

# tsne plot:
feature_df = hotsite_df.iloc[:, 6:]
target_df = hotsite_df["ideas_state"].reset_index(drop=True)

# Scaling the values:
# x = scale(feature_df.values)
# PCA analysis (PCA projection to 2D/3D):
pca_50 = PCA(n_components=50) # PCA(n_components=3)
principal_comp_50 = pca_50.fit_transform(scale(feature_df.values)) # X already normalized with fraction

time_start = time.time()
tsne = TSNE(n_components=3, verbose=1, perplexity=40, learning_rate=20, n_iter=500) # "uniq_tfcount"] >=50
tsne = TSNE(n_components=3, verbose=1, perplexity=30, learning_rate=50, n_iter=500) 
tsne_results_orig = tsne.fit_transform(feature_df.values)
# tsne_results_orig = tsne.fit_transform(principal_comp_50)
tsne_df_orig = pd.DataFrame(tsne_results_orig)

# Including 3rd dataframe for Color reference:
final_df_orig = pd.concat([tsne_df_orig, target_df], axis = 1)
final_df_orig.rename(columns={0:"tsne1", 1:"tsne2", 2:"tsne3"}, inplace=True)

# Re-map ideas states:
# final_df_orig["ideas_state"] = final_df_orig["ideas_state"].replace(replace_ideas)
map_ideas = {"Prom assoc" : "Promoter", "Strong Enh" : "Strong Enhancer", "Weak Enh" : "Weak Enhancer", "Ctcf assoc" : "CTCF assoc"}
final_df_orig["ideas_state"] = final_df_orig["ideas_state"].map(map_ideas)
values = {'ideas_state': "Other"}
final_df_orig.fillna(value=values, inplace=True)
final_df_orig.dropna(subset=["ideas_state"], inplace=True)

# Map color labels to ideas states:
map_color = {"CTCF assoc" : "purple", "Other" : "green", "Promoter" : "red", "Strong Enhancer" : "orange", "Weak Enhancer" : "yellow"}
final_df_orig["color_label"] = final_df_orig["ideas_state"].map(map_color)
final_df_orig['color_label'] = pd.Categorical(final_df_orig['color_label'])
# my_color=final_df_orig['color_label'].cat.codes

final_df_orig.to_csv(join(output_dir, "tsne_analysis_hotmotif_sites_50more_orig_wo_pca.txt"), sep="\t", header=True, index=False)
# final_df_orig["ideas_state"] = pd.Categorical(final_df_orig["ideas_state"], categories=["Weak Enhancer", "Strong Enhancer", "Promoter"], ordered=True)
print 't-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start)

plot = (ggplot(final_df_orig) + 
        aes("tsne1", "tsne2", fill="ideas_state")+ 
        geom_point(size=1, alpha=0.7) +
        # geom_text(label = final_df["ideas_state"], size=6) +
        # geom_text(label=final_df["Label"], size=8)+
        ggtitle("tSNE dimensions colored by Hotmotif site associated cis-region") +
        # theme_bw() +
        scale_x_continuous(name="x-tsne") + scale_y_continuous(name = "y-tsne" ) +
        # xlab("Principal Component 1") + ylab("Principal Component 2") +
        scale_fill_manual(name="IDEAS Anno", values=["blueviolet", "black", "red", "orange", "yellow"]) +
        guides(fill=guide_legend(title="IDEAS States")) +
        theme(
        axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot_title = element_text(size=14, face="bold"),
        legend_title = element_text(size=8, face="bold")
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
    )

ggsave(plot, join(output_dir, "tsne_analysis_hotmotif_sites_30more_orig_fig.pdf"), width=7, height=7)

###################################
# 3D scatterplot in python/seaborn/matplotlib:

import seaborn as sns
sns.set_style("white")
from mpl_toolkits.mplot3d import Axes3D

# color categories:
map_color = {"CTCF assoc" : "purple", "Other" : "green", "Promoter" : "red", "Strong Enhancer" : "orange", "Weak Enhancer" : "yellow"}
final_df_orig["color_label"] = final_df_orig["ideas_state"].map(map_color)
final_df_orig['color_label'] = pd.Categorical(final_df_orig['color_label'])
my_color=final_df_orig['color_label'].cat.codes
# final_df_orig = final_df_orig.drop('ideas_state', 1)
# my_color = ["blueviolet", "black", "red", "orange", "yellow"]

# Plot initialisation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(final_df_orig['tsne1'], final_df_orig['tsne2'], final_df_orig['tsne3'], c=my_color, cmap="Set2_r", s=10, alpha=0.6)
 
# make simple, bare axis lines through space:
xAxisLine = ((min(final_df_orig['tsne1']), max(final_df_orig['tsne1'])), (0, 0), (0,0))
ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')

yAxisLine = ((0, 0), (min(final_df_orig['tsne2']), max(final_df_orig['tsne2'])), (0,0))
ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')

zAxisLine = ((0, 0), (0,0), (min(final_df_orig['tsne3']), max(final_df_orig['tsne3'])))
ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
 
# label the axes
ax.set_xlabel("tsne1")
ax.set_ylabel("tsne2")
ax.set_zlabel("tsne3")
ax.set_title("t-SNE on Hotmotif sites")
ax.legend()

plt.savefig(join(output_dir, "t-SNE_3d_matplotlib_fig.pdf"))


##############################
# Stacked barplot with Hotmotif-sites
# Prom-Enh and other associated:

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")
output_dir = expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis")
hotmotif_file = join(output_dir, "Hotmotif_sites.bed")
loci_df = pd.read_csv(hotmotif_file, sep="\t"); loci_df.head()
select_cols = ["chrom", "chromStart", "chromEnd", "uniq_tfcount", "merged_hotmotif_id", "id"]
final_loci_df = loci_df.loc[:, select_cols]

# Assign IDEAS annotation to sorted and all tf combined fimo coords:
retain_cols = [0,1,2,3,4,5]
hotmotif_ideas_df = assign_IDEAS_State(final_loci_df, retain_cols)
hotmotif_ideas_df = hotmotif_ideas_df.sort_values(["uniq_tfcount"]).reset_index(drop=True)
# hotmotif_ideas_df1 = pd.read_pickle(join(output_dir, "Hotmotif_sites_ideas.pkl"))

# Remap ideas annotation for final plot:
final_df = hotmotif_ideas_df.copy()
final_df["ideas_state"] = final_df["ideas_state"].replace(replace_ideas)
map_ideas = {"Prom assoc" : "Promoter", "Strong Enh" : "Strong Enhancer", "Weak Enh" : "Weak Enhancer", "Ctcf assoc" : "CTCF assoc"}
final_df["ideas_state"] = final_df["ideas_state"].map(map_ideas)
values = {'ideas_state': "Other"}
final_df.fillna(value=values, inplace=True)
final_df.dropna(subset=["ideas_state"], inplace=True)

select_cols = ["uniq_tfcount", "ideas_state"]
final_anno_df = final_df.loc[:,select_cols]

""" Details for 1 or more bound sites """
final_anno_df["counter"] = 1
final_anno_df_gp = final_anno_df.groupby(["uniq_tfcount", "ideas_state"]).apply(lambda x: x["counter"].sum())
final_anno_df_reset = final_anno_df_gp.reset_index(name="element_count")
final_anno_df_reset["total"] = final_anno_df_reset.groupby(["uniq_tfcount"])["element_count"].transform(sum)
final_anno_df_reset["element_perc"] = (final_anno_df_reset["element_count"]/final_anno_df_reset["total"])*100
final_anno_df_reset["element_perc"] = final_anno_df_reset["element_perc"].astype(float).round(2)

# Order the factor/categorical variable to color legend accordingly:
final_anno_df_reset["ideas_state"] = pd.Categorical(final_anno_df_reset["ideas_state"], categories=["Promoter", "Strong Enhancer", "Weak Enhancer","CTCF assoc", "Other"], ordered=True)
map_color = {"CTCF assoc" : "purple", "Other" : "green", "Promoter" : "red", "Strong Enhancer" : "orange", "Weak Enhancer" : "yellow"}
final_anno_df_reset["color_label"] = final_anno_df_reset["ideas_state"].map(map_color)
final_anno_df_reset.to_csv(join(output_dir, "Hotmotif_sites_annotation_and_fraction.txt"), sep="\t", header=True, index=False)


#############################
# GRO-seq analysis for eRNA activity:
import pandas as pd
import pybedtools
from os.path import join, expanduser

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")
output_dir = expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis")
hotmotif_file = join(output_dir, "Hotmotif_sites.bed")
loci_df = pd.read_csv(hotmotif_file, sep="\t"); loci_df.head()
select_cols = ["chrom", "chromStart", "chromEnd", "uniq_tfcount"]
final_loci_df = loci_df.loc[:, select_cols]
loci_bedfile = pybedtools.BedTool.from_dataframe(final_loci_df)

input_dir = expanduser('~/Dropbox/for_genemodels/motifs_compinfo_analysis/Gro_seq')
gros_df = pd.read_csv(join(input_dir, "GSM2428726_HepG2_GRO-Seq.bedGraph"), sep="\t", skiprows=1, header=None)
gros_df.columns = ["chr", "start", "end", "score"]

# Getting rid of some weird line with bedtrack and others:

gros_filt_df1 = gros_df.loc[5338451:]
gros_filt_df1[["start", "end"]] = gros_filt_df1[["start", "end"]].astype(int)
gros_bedfile1 = pybedtools.BedTool.from_dataframe(gros_filt_df1)

gros_filt_df2 = gros_df.loc[:5338449]
gros_filt_df2[["start", "end"]] = gros_filt_df2[["start", "end"]].astype(int)
gros_bedfile2 = pybedtools.BedTool.from_dataframe(gros_filt_df2)


gros_df["start"].fillna("test", inplace=True)
gros_filt_df = gros_df[~(gros_df["start"] == "test")]
gros_filt_df[["start", "end"]] = gros_filt_df[["start", "end"]].astype(int)
gros_filt_sort_df = gros_filt_df.sort_values(["chr", "start", "end"])
gros_filt_sort_df = gros_filt_sort_df.reset_index(drop=True)
gros_bedfile = pybedtools.BedTool.from_dataframe(gros_filt_sort_df)

loci_gros_intersect = loci_bedfile.intersect(gros_bedfile2, wa=True, wb=True)
loci_gros_df = pd.read_csv(loci_gros_intersect.fn, sep="\t", header=None)
loci_gros_df.columns = final_loci_df.columns.tolist() + gros_df.columns.tolist()

# bin and group the scores:
bins = [1,5,10,20,30,40,50,70,100,500]
names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
loci_gros_df['binned_tf_count'] = pd.cut(loci_gros_df["uniq_tfcount"], bins, labels=names, right=False)
loci_gros_df_groupmean = loci_gros_df.groupby(["binned_tf_count"])["score"].apply(lambda x: x.mean())
loci_gros_df_groupmean1 = loci_gros_df.groupby(["uniq_tfcount"])["score"].apply(lambda x: x.mean())
print(loci_gros_df_groupmean)
print(loci_gros_df_groupmean1)


########################################
# Conservation - GERP score analysis:

import pandas as pd
import os
import re
from scipy import stats
import pybedtools
from pybedtools import BedTool
from glob import glob
from os.path import join
from os.path import splitext
from os.path import basename, expanduser

input_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/gerp_analysis/files/hg19_GERP_elements"
orig_gerp_elem_file = expanduser(join(input_dir, "hg19*"))
gerp_file_list = glob(orig_gerp_elem_file)

# Gerp element file processing:
gerp_df_list = []
for each_file in gerp_file_list:
    regex_pat = re.compile(r'hg19_(.*)_elems.txt')
    chrom = regex_pat.findall(basename(each_file))[0]
    df = pd.read_csv(each_file, sep="\t", header=None)
    df.columns = ["start", "end", "length", "score", "pval"]
    df["chrom"] = chrom
    df["strand"] = "."
    select_cols = ["chrom", "start", "end", "length", "score", "strand", "pval"]
    final_df  = df.loc[:,select_cols]
    gerp_df_list.append(final_df)

combined_gerp_df = pd.concat(gerp_df_list, ignore_index=True)
combined_gerp_sorted_df = combined_gerp_df.sort_values(["chrom","start","end"]).reset_index(drop=True)
combined_gerp_sorted_df.to_csv(join(input_dir, "combined_hg19_gerp_elems.txt"), sep="\t", header=True, index=False)
gerpelem_bedfile = pybedtools.BedTool.from_dataframe(combined_gerp_sorted_df) 

# HOTmotif site processing:
output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")
output_dir = expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis")
hotmotif_file = join(output_dir, "Hotmotif_sites.bed")
loci_df = pd.read_csv(hotmotif_file, sep="\t"); loci_df.head()
select_cols = ["chrom", "chromStart", "chromEnd", "uniq_tfcount", "binned_tf_count", "binned_tf_count_1"]
final_loci_df = loci_df.loc[:, select_cols]
loci_bedfile = pybedtools.BedTool.from_dataframe(final_loci_df)

# gerp element intersection with hotsite:
gerp_loci_intersect = loci_bedfile.intersect(gerpelem_bedfile, wao=True) 
gerp_elem_df = pd.read_csv(gerp_loci_intersect.fn, sep="\t", header=None)
# gerp_elem_df.columns = final_loci_df.columns.tolist() + final_df.columns.tolist() + ["dist_overlap"]

gerp_elem_df.rename(columns={13:"overlap_bases"}, inplace=True)
gerp_elem_df.rename(columns={10:"gerp_elem_sum_score"}, inplace=True)
gerp_elem_df.rename(columns={9:"gerp_elem_length"}, inplace=True)
gerp_elem_df.loc[gerp_elem_df["gerp_elem_length"] == ".", "gerp_elem_length"] = 0
gerp_elem_df.loc[gerp_elem_df["gerp_elem_sum_score"] == ".", "gerp_elem_sum_score"] = 0

gerp_elem_df["gerp_elem_score_perbase"] = gerp_elem_df["gerp_elem_sum_score"].astype(float)/gerp_elem_df["gerp_elem_length"].astype(float)
gerp_elem_df["gerp_elem_score_perbase"].fillna(0,inplace=True)
gerp_elem_df["overlap_score"] = gerp_elem_df["overlap_bases"] * gerp_elem_df["gerp_elem_score_perbase"]
# mask_gt = gerp_elem_df["overlap_bases"] > 0
# gerp_elem_df.loc[mask_gt, "overlap_bases"] = gerp_elem_df["overlap_bases"] # to adjust the bedtool 0 based coordinates: 

gerp_elem_df_grouped = gerp_elem_df.groupby([0,1,2,3,4])["overlap_bases", "overlap_score"].apply(lambda x: x.sum())
gerp_elem_df_overlap = gerp_elem_df_grouped.reset_index()
gerp_elem_df_overlap["total_length"] = gerp_elem_df_overlap.iloc[:,2].astype(int) - gerp_elem_df_overlap.iloc[:,1].astype(int)

gerp_elem_df_overlap["fraction_overlap_per_100bp"] = (gerp_elem_df_overlap["overlap_bases"]/gerp_elem_df_overlap["total_length"])*100
gerp_elem_df_overlap["fraction_overlap_per_kb"] = (gerp_elem_df_overlap["overlap_bases"]/gerp_elem_df_overlap["total_length"])*1000
gerp_elem_df_overlap["fraction_overlap_score_per_100bp"] = (gerp_elem_df_overlap["overlap_score"]/gerp_elem_df_overlap["total_length"])*100
gerp_elem_df_overlap["fraction_overlap_score_per_kb"] = (gerp_elem_df_overlap["overlap_score"]/gerp_elem_df_overlap["total_length"])*1000

# bins = [1,5,10,20,30,40,50,70,100,150,500]
# names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100-149", "150+" ]
# gerp_elem_df_overlap['binned_tf_count_zoomed'] = pd.cut(gerp_elem_df_overlap.iloc[:,3], bins, right=False, labels=names)

# bin based on uniq_tfcounts:
bins = [1,5,10,20,30,40,50,70,100,500]
names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
gerp_elem_df_overlap['binned_tf_count'] = pd.cut(gerp_elem_df_overlap.iloc[:,3], bins, right=False, labels=names)

gerp_elem_df_sorted = gerp_elem_df_overlap.sort_values([4]).reset_index(drop=True)
gerp_elem_df_sorted.rename(columns={0:"chrom"}, inplace=True)
gerp_elem_df_sorted.rename(columns={1:"start"}, inplace=True)
gerp_elem_df_sorted.rename(columns={2:"end"}, inplace=True)
gerp_elem_df_sorted.rename(columns={3:"uniq_tfcount"}, inplace=True)
gerp_elem_df_sorted.rename(columns={4:"misc_tfcount"}, inplace=True)
gerp_elem_df_sorted.to_csv(join(output_dir, "combined_gerp_element_analysis_boxplot_data_with_gerp_elem_overlap_fraction_final.txt"), sep="\t", header=True, index=False)



""" Perform 2 sample KS test to see if the distribution b/w 2 samples is significantly different:"""
input_file = os.path.expanduser("~/Dropbox/encode_3/gerp_analysis/analysis_output/combined_gerp_element_analysis_boxplot_data_with_gerp_elem_overlap_fraction_final_edited.txt")
gerp_df = pd.read_csv(input_file, sep="\t")

""" Boundary of seperation for low and highly bound regions """
tf_number = 80

gerp_df_lowly_bound = gerp_df.loc[gerp_df["TF_bound"] <= tf_number]
gerp_df_highly_bound = gerp_df.loc[~(gerp_df["TF_bound"] <= tf_number)]

gerp_df.shape[0] == gerp_df_lowly_bound.shape[0] + gerp_df_highly_bound.shape[0]
# gerp_lowly_bound_array = gerp_df_lowly_bound["TF_bound"].as_matrix()
# gerp_highly_bound_array = gerp_df_highly_bound["TF_bound"].as_matrix()


""" KS (Kolmogorov-Smirnov) 2 sample test for promoter associated """
gerp_df_lowly_bound_promoter = gerp_df_lowly_bound.loc[gerp_df_lowly_bound["re_anno"] == "Promoter_assoc"]
gerp_df_highly_bound_promoter = gerp_df_highly_bound.loc[gerp_df_highly_bound["re_anno"] == "Promoter_assoc"]

gerp_df_lowly_bound_promoter_array = gerp_df_lowly_bound_promoter["fraction_overlap_score_per_100bp"].as_matrix()
gerp_df_highly_bound_promoter_array = gerp_df_highly_bound_promoter["fraction_overlap_score_per_100bp"].as_matrix()

### KS test:
ks_2_sample_test  = stats.ks_2samp(gerp_df_lowly_bound_promoter_array, gerp_df_highly_bound_promoter_array)
ks_pval_prom = ks_2_sample_test[1]

if ks_pval_prom == 0:
    new_ks_pval_prom = 1e-323
    print "KS pvalue for promoter assoc :", new_ks_pval_prom
else:
    new_ks_pval_prom = ks_pval_prom
    print "KS pvalue for promoter assoc :", new_ks_pval_prom



""" KS (Kolmogorov-Smirnov) 2 sample test for enhancer associated """
gerp_df_lowly_bound_enhancer = gerp_df_lowly_bound.loc[gerp_df_lowly_bound["re_anno"] == "Enhancer_assoc"]
gerp_df_highly_bound_enhancer = gerp_df_highly_bound.loc[gerp_df_highly_bound["re_anno"] == "Enhancer_assoc"]

gerp_df_lowly_bound_enhancer_array = gerp_df_lowly_bound_enhancer["fraction_overlap_score_per_100bp"].as_matrix()
gerp_df_highly_bound_enhancer_array = gerp_df_highly_bound_enhancer["fraction_overlap_score_per_100bp"].as_matrix()

### KS test:
ks_2_sample_test  = stats.ks_2samp(gerp_df_lowly_bound_enhancer_array, gerp_df_highly_bound_enhancer_array)
ks_pval_enh = ks_2_sample_test[1]

if ks_pval_enh == 0:
    new_ks_pval_enh = 1e-323
    print "KS pvalue for enhancer assoc :", new_ks_pval_enh
else:
    new_ks_pval_enh = ks_pval_enh
    print "KS pvalue for enhancer assoc :", new_ks_pval_enh

### Final check if all the data points were included:
gerp_df_lowly_bound.shape[0] == gerp_df_lowly_bound_promoter.shape[0] + gerp_df_lowly_bound_enhancer.shape[0]
gerp_df_highly_bound.shape[0] == gerp_df_highly_bound_promoter.shape[0] + gerp_df_highly_bound_enhancer.shape[0]
gerp_df.shape[0] == gerp_df_highly_bound_promoter.shape[0] + gerp_df_lowly_bound_promoter.shape[0] + gerp_df_highly_bound_enhancer.shape[0] + gerp_df_lowly_bound_enhancer.shape[0] 


#########################################
#########################################
# Random Forest train for Ranger:

final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10_binary.pkl"))
# final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))
# final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))

final_pca_df.to_csv(join(output_dir, "Binary_Hotsite_10_TF_random_forest_train.txt"), sep="\t", header=True, index=False)

# Randomly sample 
final_pca_df_sample = final_pca_df.groupby('ideas_state').apply(lambda x: x.sample(2000, random_state=20)).reset_index(drop=True)
for key,value in final_pca_df.groupby('ideas_state'):
    print key, value.shape

# tsne plot:
feature_df = final_pca_df_sample.iloc[:, 6:]
target_df = final_pca_df_sample["ideas_state"]

# Separating out the features and convert to numpy arrays:
x = feature_df.values 

# Separating out the target or labels: 
y = target_df.values


