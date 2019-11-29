import pandas as pd, numpy as np
import seaborn as sns; sns.set()
import os, re, collections
import pybedtools
from glob import glob
import json, pickle
from os.path import splitext
from os.path import basename    
from os.path import join

main_dir = "/gpfs/gpfs1/home/schhetri/fox_meme_motifs"

output_dir = join(main_dir, "occupancy_analysis")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sub_dir = join(output_dir, "cobound_files_centrimorun_dir")
if not os.path.exists(sub_dir):
    os.makedirs(sub_dir)

# Find identifier(tf_name) for each CisBP id:
infile = open("/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        print(line)

cisbp_df = pd.Series(motif_list).str.split(expand=True)
cisbp_df.columns = ["motif", "cisbp_motif_id", "tf_name"]
cisbp_df["tf_name_final"] = cisbp_df["tf_name"].str.upper().str.replace(r"_.*", "").str.replace("\(", "").str.replace("\)", "").str.replace(r"\..*", "")

# Merge CisBP id with 208 TF names:
df_read = pd.read_csv("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/tf_list_208_tfs.txt", sep="\t", header=None)
df_read.columns = ["TF_name"]
df_read["tf_name_final"] = df_read["TF_name"].str.replace("_human|_v.*|_iso.*", "").str.upper()
df_final = df_read.sort_values(["tf_name_final"]).drop_duplicates(["tf_name_final"]).reset_index(drop=True)

merged_df = pd.merge(df_final.loc[:,["tf_name_final"]], cisbp_df.loc[:,["cisbp_motif_id", "tf_name_final"]], on=["tf_name_final"])
merged_df.to_csv("tf_list_208_tfs_withCisBP_id.txt", sep="\t", header=False, index=False)

# CisBP tf list equiv to grepping:
infile = open("/gpfs/gpfs1/home/schhetri/Tools/final_custom_CisBP_database/Homosapiens_custom_cisbp_allmotifs.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        print(line)

cisbp_df1 = pd.Series(motif_list).str.split(expand=True) # cisbp_df["a", "b", "c"] = pd.Series(motif_list).apply(pd.Series)
cisbp_df1.columns = ["motif", "cisbp_motif_id", "tf_name"]
cisbp_df1["tf_name_final"] = cisbp_df1["tf_name"].str.upper().str.replace(r"\..*", "")

infile = open("/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        print(line)

cisbp_df = pd.Series(motif_list).str.split(expand=True) # cisbp_df["a", "b", "c"] = pd.Series(motif_list).apply(pd.Series)
cisbp_df.columns = ["motif", "cisbp_motif_id", "tf_name"]
cisbp_df["tf_name_final"] = cisbp_df["tf_name"].str.upper().str.replace(r"_.*", "").str.replace("\(", "").str.replace("\)", "").str.replace(r"\..*", "")

# Find xls file with partner TFs and orthogonal TFs:
xls_file="/gpfs/gpfs1/home/schhetri/fox_meme_motifs/Forkhead_partner_TFs.xlsx"
df = pd.read_excel(xls_file, sheetname="Sheet1", header=None)
df.columns = ["primary_tf", "sec_tf", "status", "motif"]
nonfox_df = df[~df["primary_tf"].str.contains("FOX")].reset_index(drop=True)
nonfox_df["primary_tf"] = nonfox_df["primary_tf"].str.replace("NCOA2", "NCoA2") # Weird naming convention
nonfox_tf_list = nonfox_df["primary_tf"].unique().tolist()

# Find list of motif id for orthogonal factors bearing secondary fox motif:
nonfox_tf_noPwms = nonfox_df[~nonfox_df["primary_tf"].isin(cisbp_df["tf_name_final"])]
nonfox_tf_noPwms.to_csv(join(main_dir, "orthogonal_tfs_bearing_foxmotif_noPwms.txt"), sep="\t", header=True, index=False)
#nonfox_tf_withPwms = nonfox_df[nonfox_df["primary_tf"].isin(cisbp_df["tf_name_final"])]
nonfox_tf_withPwms = cisbp_df[cisbp_df["tf_name_final"].isin(nonfox_df["primary_tf"])].reset_index(drop=True)
select_cols = ["cisbp_motif_id", "tf_name", "tf_name_final"]
nonfox_tf_withPwms = nonfox_tf_withPwms.loc[:,select_cols]
nonfox_tf_withPwms.to_csv(join(main_dir, "orthogonal_tfs_bearing_foxmotif_withPwms.txt"), sep="\t", header=True, index=False)

# Generate dict for peak-bed file list :
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
file_list = [glob(join(file_pat, "*" + each + "*"))[0] for each in nonfox_tf_list]
regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_dict = {regex_pat.findall(file)[0]:file for file in file_list}
bed_dict.update({"SOX13_iso1":"/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL161284.filt.nodup.srt.SE_VS_SL161344.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_SOX13_iso1[FLAG]"})

# Generate dict for fimo-bed file list :
fimodict_dir = "/gpfs/gpfs1/home/schhetri/fox_meme_motifs/fox_peak_files/fimo_fox_peak_files"
fimopath = glob(join(fimodict_dir, "FOX*", "fimo_dir", "*.txt"))
regex_pat = re.compile("fimo_fox_peak_files\\/(.*)\\/fimo_dir")
fimo_dict = {regex_pat.findall(each)[0]:each for each in fimopath}


def combined_tf_bedfiles_fromdict(file_dict):
    concat_file_list = []
    for key,val in file_dict.iteritems():
        tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(val))[0]
        df_read = pd.read_csv(val, header=None, sep="\t")
        df_read = df_read.iloc[:,[0,1,2]]
        df_read.columns = ["chr", "start", "end"]
        df_read["tf_name"] = tf_name
        df_read["tf_peakcount"] = df_read.shape[0]
        select_cols = ["chr","start","end", "tf_name", "tf_peakcount"]
        df_read = df_read.loc[:,select_cols]
        print df_read.shape
        concat_file_list.append(df_read)

    combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
    combined_tf_df_sort = combined_tf_df.sort_values(select_cols).reset_index(drop=True)
    combined_tf_df_sort.to_csv(join(output_dir, "final_37tfs_combined_coords.bed"), header=False, index=False, sep="\t")
    return(combined_tf_df_sort)

combined_tfbed_df = combined_tf_bedfiles_fromdict(bed_dict)


def combined_tf_fimofiles_fromdict(file_dict):
    concat_file_list = []
    #file_dict = fimo_dict
    for key,value in file_dict.iteritems():
        df = pd.read_csv(value, sep="\t")
        df = df.iloc[:,1:]
        df.rename(columns={"sequence name" : "chrom"}, inplace=True)
        df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
        df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
        df["cobind_tf_name"] = key
        df["fimo_peakcount"] = df.shape[0] 
        select_cols = ["chrom", "chromStart", "chromEnd", "cobind_tf_name", "fimo_peakcount"]
        select_df = df.loc[:,select_cols]
        print("Current dimension of motif model {} : {}".format(key, select_df.shape)) 
        select_df = select_df.drop_duplicates() # select_df.duplicated(keep=False)
        print("Dropping duplicates if any, current dimension of motif model : {}\n".format(select_df.shape))
        motif_final_df = select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)
        concat_file_list.append(motif_final_df)

    combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
    combined_tf_df_sort = combined_tf_df.sort_values(select_cols).reset_index(drop=True)
    combined_tf_df_sort.to_csv(join(output_dir, "fox_combined_motif_coords.bed"), header=False, index=False, sep="\t")
    return(combined_tf_df_sort)

combined_motif_df = combined_tf_fimofiles_fromdict(fimo_dict)
            

def tf_peaks_cobind_percent_overlap(combined_tf_df_file, combined_motif_df_file):
    #combined_tf_df_file = combined_tfbed_df
    #combined_motif_df_file = combined_motif_df
    colname_df1 = combined_tf_df_file.columns.tolist()
    colname_df2 = combined_motif_df_file.columns.tolist()
    combined_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_tf_df_file)
    combined_motif_bedfile = pybedtools.BedTool.from_dataframe(combined_motif_df_file)
    tf_tf_intersect = combined_tf_bedfile.intersect(combined_motif_bedfile, wa=True, wb=True)

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["tf_name","cobind_tf_name"])

    master_dict = {} # read out
    master_offsetDist_df_list = [] # read out
    empty_offsetDist_df_dict = {}
    tf_cobind_list = []
    tf_cobind_frac = []
    for key,val in pybed_df_group:
        """ Segregate the group using groupby (with TF and co-partner) """
        peak_df = pybed_df_group.get_group(key)

        # Compute co-occupancy overlap fraction:
        unique_peak_df = peak_df.drop_duplicates(colname_df1) # Eliminate peak duplication
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, unique_peak_df.shape))
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        unique_peak_df["tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"])
        tf_cobind_list.append(key)
        tf_cobind_frac.append(tf_cobind_fraction)
        master_dict[key] = unique_peak_df

    # Combine co-occupancy overlap fraction for all factors:
    tf_cobind_df = pd.Series(tf_cobind_list).apply(pd.Series)
    tf_cobind_frac_df = pd.Series(tf_cobind_frac)
    combined_df = pd.concat([tf_cobind_df, tf_cobind_frac_df], axis=1)
    combined_df.columns = ["tf_name", "cobind_tf_name", "tf_cobind_frac"]
    heatmap_df = combined_df.pivot(index="tf_name", columns="cobind_tf_name", values="tf_cobind_frac")
    #heatmap_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif.txt"), sep="\t", header=True, index=True)
    #heatmap_df.to_excel(join(output_dir, "Factors_occupancy_with_fox_motif.xlsx"), header=True, index=True)
    sorted_df = combined_df.sort_values(["tf_cobind_frac"], ascending=False)
    #sorted_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif_sorted_meltVersion.txt"), sep="\t", header=True, index=False)

    return((heatmap_df, sorted_df))

heatmap_df_new, sorted_df_new = tf_peaks_cobind_percent_overlap(combined_tfbed_df, combined_motif_df)


def tf_peaks_cobind_filesave_forCentrimoRun(combined_tf_df_file, combined_motif_df_file, motif_extend):
    #combined_tf_df_file = combined_tfbed_df
    #combined_motif_df_file = combined_motif_df
    colname_df1 = combined_tf_df_file.columns.tolist()
    colname_df2 = combined_motif_df_file.columns.tolist()
    combined_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_tf_df_file)
    combined_motif_bedfile = pybedtools.BedTool.from_dataframe(combined_motif_df_file)
    tf_tf_intersect = combined_tf_bedfile.intersect(combined_motif_bedfile, wa=True, wb=True)

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["tf_name","cobind_tf_name"])

    master_dict = {} # read out
    for key,val in pybed_df_group:
        """ Segregate the group using groupby (with TF and co-partner) """
        peak_df = pybed_df_group.get_group(key)
        select_cols = colname_df2
        # motif_extend = 50
        cobound_peak = peak_df.loc[:,select_cols]
        cobound_peak = cobound_peak.drop_duplicates(["chrom","chromStart", "chromEnd"])
        cobound_peak["start"] = (cobound_peak["chromStart"] - motif_extend) # upstr_x2 = 40
        cobound_peak["end"] = (cobound_peak["chromEnd"] + motif_extend) + 1 # upstr_x1 =1, and + 1 for bedtool exclusivity
        filekey_name = "_".join(key).replace("_human", "").replace("[FLAG]", "")
        cobound_peak["tf_motif"] = filekey_name
        select_cols = ["chrom", "start", "end", "tf_motif"]
        cobound_peak = cobound_peak.loc[:,select_cols]
        master_dict[filekey_name] = cobound_peak
        cobound_peak.to_csv(join(sub_dir, filekey_name + "_motif_cobound.bed"), sep="\t", header=False, index=False)

    combined_df = pd.concat(master_dict).reset_index().rename(columns={"level_0":"tf_motif_name"}).drop(["level_1"], axis=1)
    select_cols = ["chrom", "start", "end", "tf_motif_name"]
    combined_df_final = combined_df.loc[:,select_cols]
    combined_df_final.to_csv(join(sub_dir, "all_factor_fox_motif_cobound_combined.txt"), sep="\t", header=True, index=False)

    return(combined_df_final)

combined_tf_foxbound_df = tf_peaks_cobind_filesave_forCentrimoRun(combined_tfbed_df, combined_motif_df, 40)


def tf_peaks_cobind_percent_overlap_withoffset_dist(combined_tf_df_file, combined_motif_df_file):
    #combined_tf_df_file = combined_tfbed_df
    #combined_motif_df_file = combined_motif_df
    colname_df1 = combined_tf_df_file.columns.tolist()
    colname_df2 = combined_motif_df_file.columns.tolist()
    combined_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_tf_df_file)
    combined_motif_bedfile = pybedtools.BedTool.from_dataframe(combined_motif_df_file)
    tf_tf_intersect = combined_tf_bedfile.intersect(combined_motif_bedfile, wa=True, wb=True)
    # print tf_tf_intersect.head()

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["tf_name","cobind_tf_name"])

    master_dict = {} # read out
    master_offsetDist_df_list = [] # read out
    empty_offsetDist_df_dict = {}
    tf_cobind_list = []
    tf_cobind_frac = []
    for key,val in pybed_df_group:
        """ Segregate the group using groupby (with TF and co-partner) """
        peak_df = pybed_df_group.get_group(key)

        # Compute co-occupancy overlap fraction:
        unique_peak_df = peak_df.drop_duplicates(colname_df1) # Eliminate peak duplication
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, unique_peak_df.shape))
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        unique_peak_df["tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"])
        tf_cobind_list.append(key)
        tf_cobind_frac.append(tf_cobind_fraction)
        master_dict[key] = unique_peak_df

        # Compute offset distance from each group of peak center:
        if not peak_df.empty:
            peak_df["mid_point"] = (peak_df["start"] + peak_df["end"])/2 
            peak_df["mid_point"] = peak_df["mid_point"].astype(int)
            peak_df["motif_mid_point"] = (peak_df["chromStart"] + peak_df["chromEnd"])/2 
            peak_df["motif_mid_point"] = peak_df["motif_mid_point"].astype(int)
            peak_df["offset_dist"] = peak_df["mid_point"] - peak_df["motif_mid_point"]
            peak_df["abs_offset_dist"] = abs(peak_df["offset_dist"])
            peak_df["motif_halfsize"] = (peak_df["chromEnd"] - peak_df["chromStart"])/2
            peak_df["motif_halfsize"] = peak_df["motif_halfsize"].astype(int)

            # Find nearest motif to the peak
            peak_df.rename(columns={"chr" : "peak_chrom"}, inplace=True)
            gpd_df = peak_df.groupby(["peak_chrom", "start", "end"])
            gpd_df_list = []
            for key_new,val_new in gpd_df:
                motif_peak_gp = gpd_df.get_group(key_new)
                motif_peak_gp = motif_peak_gp.nsmallest(1, "abs_offset_dist")
                gpd_df_list.append(motif_peak_gp)
                
            motif_peak_gp_df = pd.concat(gpd_df_list, ignore_index=True)
            motif_offset_df = motif_peak_gp_df.groupby(["tf_name", "cobind_tf_name"])["abs_offset_dist"].median().reset_index()
            motif_size_df = motif_peak_gp_df.groupby(["tf_name", "cobind_tf_name"])["motif_halfsize"].median().reset_index()
            merged_motif_df = pd.merge(motif_offset_df, motif_size_df, on=["tf_name", "cobind_tf_name"])
            merged_motif_df["adj_offset_dist"] = merged_motif_df["abs_offset_dist"] - merged_motif_df["motif_halfsize"]
            merged_motif_df.loc[merged_motif_df["adj_offset_dist"] < 0, "adj_offset_dist"] = 0
            master_offsetDist_df_list.append(merged_motif_df)
        
        else:
            print("\n{} : Motif coordinates NAN".format(key))
            empty_offsetDist_dict[key] = np.nan
            continue

    # Combine co-occupancy overlap fraction for all factors:
    tf_cobind_df = pd.Series(tf_cobind_list).apply(pd.Series)
    tf_cobind_frac_df = pd.Series(tf_cobind_frac)
    combined_df = pd.concat([tf_cobind_df, tf_cobind_frac_df], axis=1)
    combined_df.columns = ["tf_name", "cobind_tf_name", "tf_cobind_frac"]
    heatmap_df = combined_df.pivot(index="tf_name", columns="cobind_tf_name", values="tf_cobind_frac")
    heatmap_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif.txt"), sep="\t", header=True, index=True)
    heatmap_df.to_excel(join(output_dir, "Factors_occupancy_with_fox_motif.xlsx"), header=True, index=True)
    sorted_df = combined_df.sort_values(["tf_cobind_frac"], ascending=False)
    sorted_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif_sorted_meltVersion.txt"), sep="\t", header=True, index=False)

    # Combine off peak distance dataframes for all factors:
    offsetdist_df = pd.concat(master_offsetDist_df_list, ignore_index=True)
    offsetdist_heatmap_df = offsetdist_df.pivot_table(index="tf_name", columns="cobind_tf_name", values="adj_offset_dist")
    offsetdist_heatmap_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif_offsetdist.txt"), sep="\t", header=True, index=True)
    offsetdist_heatmap_df.to_excel(join(output_dir, "Factors_occupancy_with_fox_motif_offsetdist.xlsx"), header=True, index=True)
    offsetdist_sorted_df = offsetdist_df.sort_values(["adj_offset_dist"], ascending=False)
    offsetdist_sorted_df.to_csv(join(output_dir, "Factors_occupancy_with_fox_motif_offsetdist_sorted_meltVersion.txt"), sep="\t", header=True, index=False)

    return((heatmap_df,sorted_df,offsetdist_heatmap_df,offsetdist_sorted_df))

heatmap_df_new, sorted_df_new, offsetdist_df_new, offsetdist_sorted_df_new = tf_peaks_cobind_percent_overlap_withoffset_dist(combined_tfbed_df, combined_motif_df)


print "Script Run Complete!"


##########################################
# Fox Motif cobind Heatmap in R :
#########################################

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

# Output dir:
output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis/fox_meme_motifs/occupancy_analysis"

# Read file:
read_file <- fread(file.path(output_dir, "Factors_occupancy_with_fox_motif.txt"), sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
#read_df <- read_df[apply(read_df[,1:ncol(read_df)], MARGIN = 1, FUN = function(x) sd(x) != 0),]

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
# mat_data1 <- t(apply(mat_data1, 1, scale))

rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)


output_file_name <- file.path(output_dir, "Factors_cooccupancy_heatmap_with_fox_motifs.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ht1 <- Heatmap(mat_data, name="Cooccupancy Fraction", 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        col=color,
        column_title="FOX Motifs",
        row_title="Associated Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 5.5),
        column_names_gp = gpar(fontsize = 10), 
        na_col = "orange") 

ht1

dev.off()


##########################################
# Same Heatmap with box plot added though:
##########################################

output_file_name <- file.path(output_dir, "Factors_cooccupancy_heatmap_with_fox_motifs_and_boxplot.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data, axis = TRUE, ylim=c(0,max(mat_data)), border=FALSE, gp = gpar(fill = "red")))
ht1 <- Heatmap(mat_data, name="Cooccupancy Fraction",
        cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        col=color,
        column_title="FOX Motifs",
        row_title="Associated Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 4.5),
        column_names_gp = gpar(fontsize = 9), 
        na_col = "orange",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(2.8, "cm")) 

ht1

dev.off()


