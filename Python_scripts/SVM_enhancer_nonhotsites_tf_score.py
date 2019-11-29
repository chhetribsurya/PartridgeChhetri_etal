import re, os 
from os.path import join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists
# matplotlib.get_backend()
# matplotlib.use("Qt4agg")

# IDEAS segmentation file:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")

# Fasta indexing of reference genome fasta - if needed:
fasta_file = os.path.expanduser("~/for_chris/hg19-male.fa")
fasta_idx = Fasta(fasta_file)

# Extract meme E-value to filter statistically significant motifs:
script_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total"
if not os.path.exists(join(script_dir, "unique_TFs_motif_evals.txt")):
    os.system(join(script_dir, extract_meme_evalue.sh))


def assign_IDEAS_State(peaks_motif_df, cols_to_retain):
    ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
    ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")
    
    # file = join(output_dir, "final_tfcombined_fimocoords_with_fractionbind.bed")
    # peaks_motif_df = pd.read_csv(file, sep = "\t")    
    col_select = cols_to_retain
    peaks_motif_df = peaks_motif_df.iloc[:, col_select ]
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

    # motif_ideas_final_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.bed"), header=True, index=True, sep="\t")
    # motif_ideas_final_df.to_pickle(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.pkl"))
    # motif_ideas_distribution.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_motifsfraction_ideas_piechart_dist.txt"), header=True, index=True, sep="\t")

    return(motif_ideas_final_df)



#######################################################################################
## Downstream analysis using hotsites-tfbound_file cv-and-tf svmscore dictionary above:
#######################################################################################


import pandas as pd, numpy as np
import pybedtools
import os, re, pickle
from glob import glob
from os.path import join, basename

output_dir= "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"

sub_outdir= join(output_dir, "enhancer_random_nonhotsite_hotsite_factors_enrichment_5000")
if not os.path.exists(sub_outdir):
    os.makedirs(sub_outdir)

final_hotmotif_df = pd.read_pickle(join(output_dir,"Hotmotif_sites.pkl"))

 # Assign IDEAS annotation to sorted and all tf combined fimo coords:
retain_cols = [0,1,2,3,4,5,6,7,8,9]
final_hotmotif_df = assign_IDEAS_State(final_hotmotif_df, retain_cols)

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

final_hotmotif_df["ideas_anno"] = final_hotmotif_df["ideas_state"].replace(replace_ideas)
final_hotmotif_df["ideas_reanno"] = final_hotmotif_df["ideas_anno"].replace(remap_ideas)
final_hotmotif_df.rename(columns={"uniq_tfcount" : "num_tfbound", "id" : "hotsite_idx"}, inplace=True) 
final_hotmotif_df.to_csv(join(sub_outdir, "final_hotmotif_df_with_Ideas.txt"), sep="\t", index=False, header=True)
final_hotmotif_df.to_pickle(join(sub_outdir, "final_hotmotif_df_with_Ideas.pkl"))

final_hotmotif_df = pd.read_pickle(join(sub_outdir,"final_hotmotif_df_with_Ideas.pkl"))
final_hotmotif_df_enh = final_hotmotif_df[final_hotmotif_df["ideas_anno"] == "Strong Enh"]

##############################

# For centrimo enrichment run:
# Random 5000 enhancer non-hotsites
final_hotmotif_df_enh_samp = final_hotmotif_df_enh.loc[(final_hotmotif_df_enh["num_tfbound"]>=2)&(final_hotmotif_df_enh["num_tfbound"]<=10)]
final_hotmotif_df_enh_samp = final_hotmotif_df_enh_samp.sample(5000, random_state=20)
final_hotmotif_df_enh_samp.to_csv(join(sub_outdir, "final_hotmotif_df_enhancer_nonhotsites_2_10TF.txt"), sep="\t", index=False, header=False)
final_hotmotif_df_enh_samp_bed = final_hotmotif_df_enh_samp.iloc[:,0:3]
final_hotmotif_df_enh_samp_bed = final_hotmotif_df_enh_samp_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
final_hotmotif_df_enh_samp_bed.to_csv(join(sub_outdir, "NonhotsitesEnh_2_10TF.bed"), sep="\t", index=False, header=False)

# Random 5000 enhancer whole genome:
ideas_hepg2_df = pd.read_csv("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome", sep="\t")
ideas_hepg2_df.rename(columns={"name":"ideas_state"}, inplace=True)
ideas_hepg2_df["ideas_anno"] = ideas_hepg2_df["ideas_state"].replace(replace_ideas)
ideas_hepg2_df["ideas_reanno"] = ideas_hepg2_df["ideas_anno"].replace(remap_ideas)
ideas_hepg2_df = ideas_hepg2_df.iloc[:,1:]
ideas_hepg2_df.to_csv(join(sub_outdir, "ideas_hepg2_whole_genome.txt"), sep="\t", index=False, header=True)
ideas_hepg2_df.to_pickle(join(sub_outdir, "ideas_hepg2_whole_genome.pkl"))

ideas_hepg2_df_enh = ideas_hepg2_df[ideas_hepg2_df["ideas_anno"] == "Strong Enh"]
ideas_hepg2_df_enh_samp = ideas_hepg2_df_enh.sample(5000, random_state=20)
ideas_hepg2_df_enh_samp.to_csv(join(sub_outdir, "Hepg2_enhancer_whole_genome.txt"), sep="\t", index=False, header=False)
ideas_hepg2_df_enh_samp_bed = ideas_hepg2_df_enh_samp.iloc[:,0:3]
ideas_hepg2_df_enh_samp_bed = ideas_hepg2_df_enh_samp_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
ideas_hepg2_df_enh_samp_bed.to_csv(join(sub_outdir, "RandomGenomeEnh.bed"), sep="\t", index=False, header=False)

# Non-Random 2040 i.e all Hotsites: (>=50 TF) or 5000 sites(with >=35)
final_hotmotif_df_hot = final_hotmotif_df.loc[(final_hotmotif_df["num_tfbound"]>=35)]
final_hotmotif_df_hot = final_hotmotif_df_hot.sample(5000, random_state=20)
final_hotmotif_df_hot.to_csv(join(sub_outdir, "final_hotmotif_df_enhancer_hotsites_50_moreTF.txt"), sep="\t", index=False, header=False)
final_hotmotif_df_hot_bed = final_hotmotif_df_hot.iloc[:,0:3]
final_hotmotif_df_hot_bed = final_hotmotif_df_hot_bed.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
final_hotmotif_df_hot_bed.to_csv(join(sub_outdir, "HotsitesEnh_35_moreTF.bed"), sep="\t", index=False, header=False)

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

################################################

# For hotsite problem 2 (master list df) and TF_bound >=50
master_tf_svmscore_df = pd.DataFrame(master_list)
master_tf_svmscore_df.columns = ["final_hotsite_idx", "tf_bound", "score", "id", "tf_name"]
master_tf_svmscore_df = master_tf_svmscore_df.loc[(master_tf_svmscore_df["tf_bound"]>=6)&(master_tf_svmscore_df["tf_bound"]<=10)]

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
df_rank_top_bottom.to_csv(join(output_dir,"Hotmotif_sites_problem2_histogram_data_enhbound_all_24679.txt"), header=True,  index=False, sep="\t")







