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
#merged_df.to_csv("tf_list_208_tfs_withCisBP_id.txt", sep="\t", header=False, index=False)

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
nonfox_tf_noPwms_list = nonfox_tf_noPwms["primary_tf"].unique().tolist()

#nonfox_tf_withPwms = nonfox_df[nonfox_df["primary_tf"].isin(cisbp_df["tf_name_final"])]
nonfox_tf_withPwms = cisbp_df[cisbp_df["tf_name_final"].isin(nonfox_df["primary_tf"])].reset_index(drop=True)
select_cols = ["cisbp_motif_id", "tf_name", "tf_name_final"]
nonfox_tf_withPwms = nonfox_tf_withPwms.loc[:,select_cols]
nonfox_tf_withPwms_list = nonfox_tf_withPwms["tf_name_final"].unique().tolist()
#nonfox_tf_withPwms.to_csv(join(main_dir, "orthogonal_tfs_bearing_foxmotif_withPwms.txt"), sep="\t", header=True, index=False)

# Generate nonfox tf_bed peak dict:
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
file_list = [glob(join(file_pat, "*" + each + "*"))[0] for each in nonfox_tf_list ]
regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_dict = {regex_pat.findall(file)[0]:file for file in file_list}
bed_dict.update({"SOX13_iso1":"/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL161284.filt.nodup.srt.SE_VS_SL161344.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_SOX13_iso1[FLAG]"})

# Generate nonfox tf_bed peak dict (for 7 tfs without pwms):
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
file_list = [glob(join(file_pat, "*" + each + "*"))[0] for each in nonfox_tf_noPwms_list]
regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_dict_wo_pwms = {regex_pat.findall(file)[0]:file for file in file_list}

# Generate nonfox tf_bed peak dict (for 30 tfs with pwms):
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
file_list = [glob(join(file_pat, "*" + each + "*"))[0] for each in nonfox_tf_withPwms_list]
regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_dict_with_pwms = {regex_pat.findall(file)[0]:file for file in file_list}
bed_dict_with_pwms.update({"SOX13_iso1[FLAG]":"/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL161284.filt.nodup.srt.SE_VS_SL161344.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak_SOX13_iso1[FLAG]"})

# Generate FOX peak dict:
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
file_list = glob(join(file_pat, "*" + "FOX" + "*"))
regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
fox_bed_dict = {regex_pat.findall(file)[0]:file for file in file_list}

# Generate dict for fox fimo-bed :
fimodict_dir = "/gpfs/gpfs1/home/schhetri/fox_meme_motifs/fox_peak_files/fimo_fox_peak_files"
fimopath = glob(join(fimodict_dir, "FOX*", "fimo_dir", "*.txt"))
regex_pat = re.compile("fimo_fox_peak_files\\/(.*)\\/fimo_dir")
fox_fimo_dict = {regex_pat.findall(each)[0]:each for each in fimopath}

# Generate dict for nonfox tfs fimo-bed :
fimodict_dir = "/gpfs/gpfs1/home/schhetri/fox_meme_motifs/nonfox_orthogonalTF_peak_files/fimo_nonfox_peak_files"
fimopath = glob(join(fimodict_dir, "*", "fimo_dir", "fimo.txt"))
regex_pat = re.compile("fimo_nonfox_peak_files\\/(.*)\\/fimo_dir")
nonfox_fimo_dict = {regex_pat.findall(each)[0]:each for each in fimopath}


# Peak level readout:
def combined_tf_bedfiles_fromdict(file_dict):
    concat_file_list = []
    for key,val in file_dict.iteritems():
        tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(val))[0]
        df_read = pd.read_csv(val, header=None, sep="\t")
        df_read = df_read.iloc[:,[0,1,2]]
        df_read.columns = ["chr", "start", "end"]
        df_read["tf_name"] = tf_name
        df_read["tf_name_final"] = df_read["tf_name"].str.replace(r"SOX13_iso1\[FLAG\]", "SOX13.iso1").str.replace(r"\[FLAG\]|_human|_v.*|_iso.*", "") \
                                    .str.upper().str.replace(".ISO1", ".iso1")
        df_read["tf_peakcount"] = df_read.shape[0]
        select_cols = ["chr","start","end", "tf_name", "tf_name_final", "tf_peakcount"]
        df_read = df_read.loc[:,select_cols]
        print df_read.shape
        concat_file_list.append(df_read)

    combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
    combined_tf_df_sort = combined_tf_df.sort_values(select_cols).reset_index(drop=True)
    combined_tf_df_sort.to_csv(join(output_dir, "final_37tfs_combined_coords.bed"), header=False, index=False, sep="\t")
    return(combined_tf_df_sort)

combined_nonfox_tfbed_df = combined_tf_bedfiles_fromdict(bed_dict)
combined_fox_tfbed_df = combined_tf_bedfiles_fromdict(fox_bed_dict)

# further breaking the category with and without pwms:
combined_nonfox_tfbed_nopwms_df = combined_tf_bedfiles_fromdict(bed_dict_wo_pwms)
combined_nonfox_tfbed_withpwms_df = combined_tf_bedfiles_fromdict(bed_dict_with_pwms)


# Motif level readout :
def combined_tf_fimofiles_fromdict(file_dict):
    concat_file_list = []
    #file_dict = nonfox_fimo_dict
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

combined_fox_motif_df = combined_tf_fimofiles_fromdict(fox_fimo_dict)
combined_fox_motif_df["cobind_sym_tf_name"] = "FOX"            
combined_nonfox_self_motif_df = combined_tf_fimofiles_fromdict(nonfox_fimo_dict)

# Find the list of TFs with self motif for cobind fraction compute:
tf_list_with_selfmotif = combined_nonfox_self_motif_df["cobind_tf_name"].unique()

def tf_motif_cobind_percent_overlap(tf_list_of_int, combined_nonfox_tf_df_file, combined_nonfox_self_motif_df_file, combined_fox_motif_df_file ):
    #combined_nonfox_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_nonfox_tf_df_file)
    combined_fox_motif_bedfile = pybedtools.BedTool.from_dataframe(combined_fox_motif_df_file)

    # Grouping the motifs:
    self_motif_groupby = combined_nonfox_self_motif_df_file.groupby(["cobind_tf_name"])
    nonfox_tf_bed_groupby = combined_nonfox_tf_df_file.groupby(["tf_name_final"])

    master_dict = {}
    both_self_fox_motif_peak = []
    fox_motif_only_peak = []
    self_motif_only_peak = []
    no_motif_peak = []

    #for each in combined_nonfox_self_motif_df_file["cobind_tf_name"].unique():
    for each in tf_list_of_int:
        print("Processing {}\n".format(each))
        self_motif_df_file = self_motif_groupby.get_group(each)
        nonfox_tf_df_file = nonfox_tf_bed_groupby.get_group(each)
 
        self_motif_bedfile = pybedtools.BedTool.from_dataframe(self_motif_df_file)
        nonfox_tf_bedfile = pybedtools.BedTool.from_dataframe(nonfox_tf_df_file)

        ## First case of intersection (peak&selfmotif+ peak&Foxmotif+ i.e Both)
        tf_selfmotif_int = nonfox_tf_bedfile.intersect(self_motif_bedfile, u=True)
        tf_self_fox_motif_int = tf_selfmotif_int.intersect(combined_fox_motif_bedfile, u=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_self_fox_motif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        both_self_fox_motif_peak.append((each, tf_cobind_fraction))
        master_dict["both_self_fox_motif"] = {}
        master_dict["both_self_fox_motif"].update({each:unique_peak_df})

        ## Second case of intersection (peak&selfmotif- peak&Foxmotif+ i.e Foxmotif only)
        tf_selfmotif_int = nonfox_tf_bedfile.intersect(self_motif_bedfile, v=True)
        tf_foxmotif_int = tf_selfmotif_int.intersect(combined_fox_motif_bedfile, u=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_foxmotif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        fox_motif_only_peak.append((each, tf_cobind_fraction))
        master_dict["fox_motif_only"] = {}
        master_dict["fox_motif_only"].update({each:unique_peak_df})

        ## Third case of intersection (peak&selfmotif+ peak&Foxmotif- i.e Selfmotif only)
        tf_selfmotif_int = nonfox_tf_bedfile.intersect(self_motif_bedfile, u=True)
        tf_selfmotif_int = tf_selfmotif_int.intersect(combined_fox_motif_bedfile, v=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_selfmotif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        self_motif_only_peak.append((each, tf_cobind_fraction))
        master_dict["self_motif_only"] = {}
        master_dict["self_motif_only"].update({each:unique_peak_df})

        ## Fourth case of intersection (peak&selfmotif- peak&Foxmotif- i.e No motifs at all)
        tf_selfmotif_int = nonfox_tf_bedfile.intersect(self_motif_bedfile, v=True)
        tf_without_motif_int = tf_selfmotif_int.intersect(combined_fox_motif_bedfile, v=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_without_motif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        no_motif_peak.append((each, tf_cobind_fraction))
        master_dict["no_motif"] = {}
        master_dict["no_motif"].update({each:unique_peak_df})

    ## Combine co-occupancy overlap fraction for all factors:
    both_motif_cobind_df = pd.Series(both_self_fox_motif_peak).apply(pd.Series)
    both_motif_cobind_df.columns = ["tf_name", "Both motif"]

    foxonly_cobind_df = pd.Series(fox_motif_only_peak).apply(pd.Series)
    foxonly_cobind_df.columns = ["tf_name", "FOX only"]

    selfonly_cobind_df = pd.Series(self_motif_only_peak).apply(pd.Series)
    selfonly_cobind_df.columns = ["tf_name", "Self only"]

    nomotif_cobind_df = pd.Series(no_motif_peak).apply(pd.Series)
    nomotif_cobind_df.columns = ["tf_name", "No motif"]

    df_list = [both_motif_cobind_df, foxonly_cobind_df, selfonly_cobind_df, nomotif_cobind_df]
    combined_cobind_df = reduce( lambda left, right : pd.merge(left, right, on=["tf_name"]), df_list)
    return(combined_cobind_df)


heatmap_df_withmotif = tf_motif_cobind_percent_overlap(tf_list_with_selfmotif, combined_nonfox_tfbed_df, combined_nonfox_self_motif_df, combined_fox_motif_df)
heatmap_df_withmotif.to_excel(join(output_dir, "Factors_occupancy_with_both_fox_self_and_none_motif.xlsx"), header=True, index=False)
heatmap_df_withmotif.to_csv(join(output_dir, "Factors_occupancy_with_both_fox_self_and_none_motif.txt"), sep="\t", header=True, index=False)


# Find the list of TFs without self motif for cobind fraction compute:
all_tf_name_df = pd.Series(combined_nonfox_tfbed_df["tf_name_final"].unique())
tf_list_without_selfmotif = all_tf_name_df[~all_tf_name_df.isin(combined_nonfox_self_motif_df["cobind_tf_name"])].unique()

def tf_without_selfmotif_cobind_percent_overlap(tf_list_of_int, combined_nonfox_tf_df_file, combined_fox_motif_df_file ):
    #combined_nonfox_tf_bedfile = pybedtools.BedTool.from_dataframe(combined_nonfox_tf_df_file)
    combined_fox_motif_bedfile = pybedtools.BedTool.from_dataframe(combined_fox_motif_df_file)

    # Grouping the motifs:
    nonfox_tf_bed_groupby = combined_nonfox_tf_df_file.groupby(["tf_name_final"])

    master_dict = {}
    both_self_fox_motif_peak = []
    fox_motif_only_peak = []
    self_motif_only_peak = []
    no_motif_peak = []

    #for each in combined_nonfox_self_motif_df_file["cobind_tf_name"].unique():
    for each in tf_list_of_int:
        print("Processing {}\n".format(each))
        nonfox_tf_df_file = nonfox_tf_bed_groupby.get_group(each) 
        nonfox_tf_bedfile = pybedtools.BedTool.from_dataframe(nonfox_tf_df_file)

        ## Second case of intersection (peak&selfmotif- peak&Foxmotif+ i.e Foxmotif only)
        tf_foxmotif_int = nonfox_tf_bedfile.intersect(combined_fox_motif_bedfile, u=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_foxmotif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        fox_motif_only_peak.append((each, tf_cobind_fraction))
        master_dict["fox_motif_only"] = {}
        master_dict["fox_motif_only"].update({each:unique_peak_df})

        ## Fourth case of intersection (peak&selfmotif- peak&Foxmotif- i.e No motifs at all)
        tf_without_motif_int = nonfox_tf_bedfile.intersect(combined_fox_motif_bedfile, v=True)

        """ read intersect file """
        unique_peak_df = pd.read_csv(tf_without_motif_int.fn, sep="\t", header=None)
        unique_peak_df.columns = combined_nonfox_tf_df_file.columns.tolist()
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        no_motif_peak.append((each, tf_cobind_fraction))
        master_dict["no_motif"] = {}
        master_dict["no_motif"].update({each:unique_peak_df})

    ## Combine co-occupancy overlap fraction for all factors:
    foxonly_cobind_df = pd.Series(fox_motif_only_peak).apply(pd.Series)
    foxonly_cobind_df.columns = ["tf_name", "FOX only"]

    nomotif_cobind_df = pd.Series(no_motif_peak).apply(pd.Series)
    nomotif_cobind_df.columns = ["tf_name", "No motif"]

    df_list = [foxonly_cobind_df, nomotif_cobind_df]
    combined_cobind_df = reduce( lambda left, right : pd.merge(left, right, on=["tf_name"]), df_list)
    return(combined_cobind_df)

heatmap_df_nomotif = tf_without_selfmotif_cobind_percent_overlap(tf_list_without_selfmotif, combined_nonfox_tfbed_df, combined_fox_motif_df)
heatmap_df_nomotif["Both motif"] = np.nan
heatmap_df_nomotif["Self only"] = np.nan

# to sync the order with previous category for concat:
select_cols = ["tf_name", "Both motif", "FOX only", "Self only", "No motif"]
heatmap_df_nomotif_new = heatmap_df_nomotif.loc[:,select_cols]

merged_heatmap_df = pd.concat([heatmap_df_withmotif, heatmap_df_nomotif_new], ignore_index=True)
merged_heatmap_df = merged_heatmap_df.loc[merged_heatmap_df["tf_name"] != "SOX13"]
merged_heatmap_df = merged_heatmap_df.set_index("tf_name").reset_index()
merged_heatmap_df.to_excel(join(output_dir, "Factors_occupancy_with_motif_details_37factors_fin.xlsx"), header=True, index=False)
merged_heatmap_df.to_csv(join(output_dir, "Factors_occupancy_with_motif_details_37factors_fin.txt"), sep="\t", header=True, index=False)


def foxome_cobind_percent_overlap(combined_fox_tfbed_df_file, combined_nonfox_tfbed_df_file):
    # For distinct col name - required for dropping duplicates:
    combined_fox_tfbed_df_file = combined_fox_tfbed_df_file.add_prefix("fox_")
    # combined_fox_tfbed_df_file["fox_sym_tf_name"] = "FOX"  # if needed for cumulative fox-tfcobind fraction
    colname_df1 = combined_fox_tfbed_df_file.columns.tolist()
    colname_df2 = combined_nonfox_tfbed_df_file.columns.tolist()
    combined_nonfox_tfbed_df_bedfile = pybedtools.BedTool.from_dataframe(combined_nonfox_tfbed_df_file)
    combined_fox_motif_df_bedfile = pybedtools.BedTool.from_dataframe(combined_fox_tfbed_df_file)
    tf_tf_intersect = combined_fox_motif_df_bedfile.intersect(combined_nonfox_tfbed_df_bedfile, wa=True, wb=True)

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["fox_tf_name_final","tf_name_final"])

    master_dict = {} # read out
    tf_cobind_list = []
    tf_cobind_frac = []
    for key,val in pybed_df_group:
        """ Segregate the group using groupby(with foxTF and co-partner) """
        peak_df = pybed_df_group.get_group(key)

        # Compute co-occupancy overlap fraction:
        unique_peak_df = peak_df.drop_duplicates(colname_df1) # Eliminate peak duplication
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, unique_peak_df.shape))
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"].unique().astype(float))[0]
        unique_peak_df["tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"])
        tf_cobind_list.append(key)
        tf_cobind_frac.append(tf_cobind_fraction)
        master_dict[key] = unique_peak_df

    # Combine co-occupancy overlap fraction for all factors:
    tf_cobind_df = pd.Series(tf_cobind_list).apply(pd.Series)
    tf_cobind_frac_df = pd.Series(tf_cobind_frac)
    combined_df = pd.concat([tf_cobind_df, tf_cobind_frac_df], axis=1)
    combined_df.columns = ["tf_name", "cobind_tf_name", "tf_cobind_frac"]
    combined_df = combined_df.loc[combined_df["cobind_tf_name"] != "SOX13.iso1"]
    heatmap_df = combined_df.pivot(index="tf_name", columns="cobind_tf_name", values="tf_cobind_frac")
    sorted_df = combined_df.sort_values(["tf_cobind_frac"], ascending=False).reset_index(drop=True)
    return((heatmap_df, sorted_df))

heatmap_df_foxome, sorted_cobind_df = foxome_cobind_percent_overlap(combined_fox_tfbed_df, combined_nonfox_tfbed_df)
heatmap_df_foxome.to_csv(join(output_dir, "FOXome_occupancy_with_follow_factors.txt"), sep="\t", header=True, index=True)
heatmap_df_foxome.to_excel(join(output_dir, "FOXome_occupancy_with_follow_factors.xlsx"), header=True, index=True)

# Tag heatmap dataframe with peakcount for each factor:
heatmap_df_foxome_transp = heatmap_df_foxome.transpose().reset_index()
heatmap_df_foxome_transp.rename(columns={"cobind_tf_name": "tf_name_final"}, inplace=True)
extract_peak_count_df = combined_nonfox_tfbed_df.loc[:,["tf_name_final", "tf_peakcount"]]
extract_peak_count_df = extract_peak_count_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)

merged_heatmap_df = pd.merge(heatmap_df_foxome_transp, extract_peak_count_df, on=["tf_name_final"])
merged_heatmap_df_final = merged_heatmap_df.set_index("tf_name_final").sort_values(["tf_peakcount"]).transpose()
merged_heatmap_df_final.to_csv(join(output_dir, "FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=True, index=True)
merged_heatmap_df_final.to_excel(join(output_dir, "FOXome_occupancy_with_follow_factors_peakOrdered.xlsx"), header=True, index=True)

# Eliminate last columns with peak count info:
merged_heatmap_final = merged_heatmap_df_final.iloc[:-1]


##################################################
# Normalizing factor usage for cobinding overlap:
##################################################


def normalized_foxome_cobind_percent_overlap(combined_fox_tfbed_df_file, combined_nonfox_tfbed_df_file):
    # combined_fox_tfbed_df_file = combined_fox_tfbed_df
    # combined_nonfox_tfbed_df_file = combined_nonfox_tfbed_df

    # For distinct col name - required for dropping duplicates:
    combined_fox_tfbed_df_file = combined_fox_tfbed_df_file.add_prefix("fox_")
    # combined_fox_tfbed_df_file["fox_sym_tf_name"] = "FOX"  # if needed for cumulative fox-tfcobind fraction
    colname_df1 = combined_fox_tfbed_df_file.columns.tolist()
    colname_df2 = combined_nonfox_tfbed_df_file.columns.tolist()
    combined_nonfox_tfbed_df_bedfile = pybedtools.BedTool.from_dataframe(combined_nonfox_tfbed_df_file)
    combined_fox_motif_df_bedfile = pybedtools.BedTool.from_dataframe(combined_fox_tfbed_df_file)
    tf_tf_intersect = combined_fox_motif_df_bedfile.intersect(combined_nonfox_tfbed_df_bedfile, wa=True, wb=True)
    norm_tf_tf_intersect = combined_nonfox_tfbed_df_bedfile.intersect(combined_fox_motif_df_bedfile, wa=True, wb=True)

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["fox_tf_name_final","tf_name_final"])

    # Normalized Read intersection:
    norm_pybed_df = pd.read_csv(norm_tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df2 + colname_df1  # need to reverse as pybed intersect got flipped
    norm_pybed_df.columns = header
    norm_pybed_df_group = norm_pybed_df.groupby(["tf_name_final", "fox_tf_name_final"])

    master_dict = {} # read out
    tf_cobind_list = []
    tf_cobind_frac = []
    for key,val in pybed_df_group:
        """ Segregate the group using groupby(with foxTF and co-partner) """
        peak_df = pybed_df_group.get_group(key)
        norm_peak_df = norm_pybed_df_group.get_group(key[::-1]) # reverse key as pybed intersect got flipped

        # Compute co-occupancy overlap fraction:
        unique_peak_df = peak_df.drop_duplicates(colname_df1) # Eliminate peak duplication
        norm_unique_peak_df = norm_peak_df.drop_duplicates(colname_df2) # Eliminate peak duplication
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, unique_peak_df.shape))
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, norm_unique_peak_df.shape))
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"].unique().astype(float))[0]
        normalizing_factor = (norm_unique_peak_df.shape[0]/norm_unique_peak_df["tf_peakcount"].unique().astype(float))[0]
        norm_tf_cobind_fraction = tf_cobind_fraction * normalizing_factor
        unique_peak_df["tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"])
        unique_peak_df["norm_tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"]) * (norm_unique_peak_df.shape[0]/norm_unique_peak_df["tf_peakcount"])
        tf_cobind_list.append(key)
        tf_cobind_frac.append(norm_tf_cobind_fraction) # normalized fraction
        # tf_cobind_frac.append(tf_cobind_fraction) # normalized fraction
        master_dict[key] = unique_peak_df

    # Combine co-occupancy overlap fraction for all factors:
    tf_cobind_df = pd.Series(tf_cobind_list).apply(pd.Series)
    tf_cobind_frac_df = pd.Series(tf_cobind_frac)
    combined_df = pd.concat([tf_cobind_df, tf_cobind_frac_df], axis=1)
    combined_df.columns = ["tf_name", "cobind_tf_name", "tf_cobind_frac"]
    combined_df = combined_df.loc[combined_df["cobind_tf_name"] != "SOX13.iso1"]
    heatmap_df = combined_df.pivot(index="tf_name", columns="cobind_tf_name", values="tf_cobind_frac")
    sorted_df = combined_df.sort_values(["tf_cobind_frac"], ascending=False).reset_index(drop=True)
    return((heatmap_df, sorted_df))

norm_heatmap_df_foxome, norm_sorted_cobind_df = normalized_foxome_cobind_percent_overlap(combined_fox_tfbed_df, combined_nonfox_tfbed_df)
norm_heatmap_df_foxome.to_csv(join(output_dir, "norm_FOXome_occupancy_with_follow_factors.txt"), sep="\t", header=True, index=True)
norm_heatmap_df_foxome.to_excel(join(output_dir, "norm_FOXome_occupancy_with_follow_factors.xlsx"), header=True, index=True)

# Tag heatmap dataframe with peakcount for each factor:
norm_heatmap_df_foxome_transp = norm_heatmap_df_foxome.transpose().reset_index()
norm_heatmap_df_foxome_transp.rename(columns={"cobind_tf_name": "tf_name_final"}, inplace=True)
extract_peak_count_df = combined_nonfox_tfbed_df.loc[:,["tf_name_final", "tf_peakcount"]]
extract_peak_count_df = extract_peak_count_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)

merged_heatmap_df = pd.merge(norm_heatmap_df_foxome_transp, extract_peak_count_df, on=["tf_name_final"])
merged_heatmap_df_final = merged_heatmap_df.set_index("tf_name_final").sort_values(["tf_peakcount"]).transpose()
merged_heatmap_df_final.to_csv(join(output_dir, "norm_FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=True, index=True)
merged_heatmap_df_final.to_excel(join(output_dir, "norm_FOXome_occupancy_with_follow_factors_peakOrdered.xlsx"), header=True, index=True)


##################################################
# Normalizing factor usage for cobinding overlap:
##################################################


def mean_normalized_foxome_cobind_percent_overlap(combined_fox_tfbed_df_file, combined_nonfox_tfbed_df_file):
    # combined_fox_tfbed_df_file = combined_fox_tfbed_df
    # combined_nonfox_tfbed_df_file = combined_nonfox_tfbed_df

    # For distinct col name - required for dropping duplicates:
    combined_fox_tfbed_df_file = combined_fox_tfbed_df_file.add_prefix("fox_")
    # combined_fox_tfbed_df_file["fox_sym_tf_name"] = "FOX"  # if needed for cumulative fox-tfcobind fraction
    colname_df1 = combined_fox_tfbed_df_file.columns.tolist()
    colname_df2 = combined_nonfox_tfbed_df_file.columns.tolist()
    combined_nonfox_tfbed_df_bedfile = pybedtools.BedTool.from_dataframe(combined_nonfox_tfbed_df_file)
    combined_fox_motif_df_bedfile = pybedtools.BedTool.from_dataframe(combined_fox_tfbed_df_file)
    tf_tf_intersect = combined_fox_motif_df_bedfile.intersect(combined_nonfox_tfbed_df_bedfile, wa=True, wb=True)

    # Read intersection:
    pybed_df = pd.read_csv(tf_tf_intersect.fn, sep="\t", header=None)
    header = colname_df1 + colname_df2
    pybed_df.columns = header
    pybed_df_group = pybed_df.groupby(["fox_tf_name_final","tf_name_final"])

    # Normalize using mean peak counts for factors:
    peak_mean_compute_df = pybed_df.loc[:,["tf_name_final", "tf_peakcount"]]
    peak_mean_compute_df = peak_mean_compute_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)
    mean_peak_count = peak_mean_compute_df["tf_peakcount"].mean()

    master_dict = {} # read out
    tf_cobind_list = []
    tf_cobind_frac = []
    for key,val in pybed_df_group:
        """ Segregate the group using groupby(with foxTF and co-partner) """
        peak_df = pybed_df_group.get_group(key)

        # Compute co-occupancy overlap fraction:
        unique_peak_df = peak_df.drop_duplicates(colname_df1) # Eliminate peak duplication
        print("\nProcessing {} : initial shape {} and later shape {}".format(key, val.shape, unique_peak_df.shape))
        tf_cobind_fraction = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"].unique().astype(float))[0]
        mean_normalizing_factor = (mean_peak_count/unique_peak_df["tf_peakcount"].unique().astype(float)[0]) # mean normalization with respect to total peaks
        norm_tf_cobind_fraction = tf_cobind_fraction * mean_normalizing_factor
        unique_peak_df["tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"])
        unique_peak_df["norm_tf_cobind_frac"] = (unique_peak_df.shape[0]/unique_peak_df["fox_tf_peakcount"]) * mean_normalizing_factor
        tf_cobind_list.append(key)
        tf_cobind_frac.append(norm_tf_cobind_fraction) # normalized fraction
        # tf_cobind_frac.append(tf_cobind_fraction) # normalized fraction
        master_dict[key] = unique_peak_df

    # Combine co-occupancy overlap fraction for all factors:
    tf_cobind_df = pd.Series(tf_cobind_list).apply(pd.Series)
    tf_cobind_frac_df = pd.Series(tf_cobind_frac)
    combined_df = pd.concat([tf_cobind_df, tf_cobind_frac_df], axis=1)
    combined_df.columns = ["tf_name", "cobind_tf_name", "tf_cobind_frac"]
    combined_df = combined_df.loc[combined_df["cobind_tf_name"] != "SOX13.iso1"]
    heatmap_df = combined_df.pivot(index="tf_name", columns="cobind_tf_name", values="tf_cobind_frac")
    sorted_df = combined_df.sort_values(["tf_cobind_frac"], ascending=False).reset_index(drop=True)
    return((heatmap_df, sorted_df))

mean_norm_heatmap_df_foxome, mean_norm_sorted_cobind_df = mean_normalized_foxome_cobind_percent_overlap(combined_fox_tfbed_df, combined_nonfox_tfbed_df)
mean_norm_heatmap_df_foxome.to_csv(join(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors.txt"), sep="\t", header=True, index=True)
mean_norm_heatmap_df_foxome.to_excel(join(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors.xlsx"), header=True, index=True)

# Tag heatmap dataframe with peakcount for each factor:
mean_norm_heatmap_df_foxome_transp = mean_norm_heatmap_df_foxome.transpose().reset_index()
mean_norm_heatmap_df_foxome_transp.rename(columns={"cobind_tf_name": "tf_name_final"}, inplace=True)
extract_peak_count_df = combined_nonfox_tfbed_df.loc[:,["tf_name_final", "tf_peakcount"]]
extract_peak_count_df = extract_peak_count_df.drop_duplicates(["tf_name_final"]).reset_index(drop=True)

merged_heatmap_df = pd.merge(mean_norm_heatmap_df_foxome_transp, extract_peak_count_df, on=["tf_name_final"])
merged_heatmap_df_final = merged_heatmap_df.set_index("tf_name_final").sort_values(["tf_peakcount"]).transpose()
merged_heatmap_df_final.to_csv(join(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=True, index=True)
merged_heatmap_df_final.to_excel(join(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors_peakOrdered.xlsx"), header=True, index=True)

