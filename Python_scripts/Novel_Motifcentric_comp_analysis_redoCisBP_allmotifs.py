import re, os 
from os.path import join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists

import matplotlib
matplotlib.get_backend()
matplotlib.use("Qt4agg")
import matplotlib.pyplot as plt
from plotnine import *
#from ggplot import *

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis_tomtom_redo"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# IDEAS segmentation file:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")

# Fasta indexing of reference genome fasta - if needed:
fasta_file = os.path.expanduser("~/for_chris/hg19-male.fa")
fasta_idx = Fasta(fasta_file)

# Load TF annoation file:
# anno_df = pd.read_csv("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df = pd.read_csv("~/for_chris/batch_I/TFs_Annotation_file.txt", sep="\t")
dbf_anno_df = anno_df.loc[anno_df["Category"] == "DBF"]

# Extract meme E-value to filter statistically significant motifs:
script_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total"
if not os.path.exists(join(script_dir, "unique_TFs_motif_evals.txt")):
    os.system(join(script_dir, extract_meme_evalue.sh))

# Read motif E-val containing motif file - output of above extract_meme_evalue.sh (script)
df_motif_table = pd.read_csv(join(script_dir, "unique_TFs_motif_evals.txt"), sep="\t")
df_motif_table["E-value_sig"] = np.where(df_motif_table["E-value"].astype(float) < 0.05, df_motif_table["E-value"], "NA")

# Generate dict for peak-bed file list :
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"
bed_filelist = glob(file_pat)
bed_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_filedict = {bed_regex_pat.findall(file)[0]:file for file in bed_filelist}
# df = pd.DataFrame(bed_filedict.items())

# Generate dict for meme file list :
file_pat = join(script_dir, "unique_TFs/*_narrowPeak_*/meme_chip/meme_out/meme.txt")
meme_filelist = glob(file_pat)
meme_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*?)/meme')
meme_filedict = {meme_regex_pat.findall(file)[0]:file for file in meme_filelist}
# df1 = pd.DataFrame(meme_filedict.items())

# Generate dict for fimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_filelist = glob(file_pat)
fimo_regex_pat = re.compile(r'.*unique_TFs/.*_fimo_motif_(.*).txt')
fimo_filedict = {fimo_regex_pat.findall(file)[0]:file for file in fimo_filelist}
# df2 = pd.DataFrame(fimo_filedict.items())

# Generate dict for centrimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/unique_TFs/*_narrowPeak_*/meme_cetrimo_dir/centrimo.txt"
centrimo_filelist = glob(file_pat)
centrimo_regex_pat = re.compile(r'unique_TFs/.*_narrowPeak_(.*?)/meme_*') 
centrimo_filedict = {centrimo_regex_pat.findall(file)[0]:file for file in centrimo_filelist}
# df3 = pd.DataFrame(centrimo_filedict.items())

# Generate dict for tomtom file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/unique_TFs/*_narrowPeak_*/meme_tomtom_dir/tomtom.txt"
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/*_narrowPeak_*/meme_jaspar_tomtom_dir/tomtom.txt"
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/*_narrowPeak_*/meme_cisbp_tomtom_dir/tomtom.txt"
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/*_narrowPeak_*/meme_finalredo_cisbp_tomtom_dir/tomtom.txt"
tomtom_filelist = glob(file_pat)
tomtom_regex_pat = re.compile(r'finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/.*_narrowPeak_(.*?)/meme_*')
tomtom_filedict = {tomtom_regex_pat.findall(file)[0]:file for file in tomtom_filelist}
# df4 = pd.DataFrame(tomtom_filedict.items())


def final_fimo_motifs_model(fimo_file, tf_name, **kwargs):
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


def final_centrimo_motifs_model(centrimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    # with open(centrimo_file,"r") as file:
    # centrimo_dict = {}
    centrimo_keylist = []
    centrimo_valuelist = []
    for line in open(centrimo_file,"r").readlines(): #[:10]
        if not line.startswith("#"):
            splitted = line.strip("\n").split("\t")
            centrimo_motif_id, e_value, bin_location, bin_width = splitted[1], splitted[2], splitted[5], splitted[6]
            regex_pat = re.compile(r' +') # replace more than one space
            motif_id, alt_name = regex_pat.sub("|", centrimo_motif_id).split("|")
            centrimo_values = [alt_name, e_value, bin_location, bin_width]
            centrimo_keylist.append(motif_id)
            centrimo_valuelist.append(centrimo_values) #centrimo_dict[motif_id] = centrimo_values
    

    df = pd.concat([pd.Series(centrimo_keylist), pd.Series(centrimo_valuelist)], axis=1)
    if not df.empty:
        df.columns = ["motif_id", "altname"]
        new_col_names = ["altname", "E-value", "bin_location", "bin_width"] 
        df[new_col_names] = df["altname"].apply(pd.Series)

        # Setting values in dataframe using .loc:
        motif_idlist = map(str, motif_list)
        df.loc[df["motif_id"].isin(motif_idlist), "motif_id"] = "MOTIF" + df["motif_id"]
        
        # Rename motifs id to select motifs requested in args:
        select_motif_list = ["MOTIF" + each for each in motif_idlist]
        final_df = df.loc[df["motif_id"].isin(select_motif_list)]

    else:
        new_col_names = ["altname", "E-value", "bin_location", "bin_width"] 
        final_df = pd.DataFrame({"motif_id":["MOTIF1"], "altname":[None], "E-value":[None], "bin_location":[None], "bin_width":[None]})
        select_cols = ["motif_id"] + new_col_names # order cols
        final_df = final_df.loc[:, select_cols]

    return(final_df)


def final_tomtom_motifs_model(tomtom_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(tomtom_file, sep="\t")
    df.rename(columns={"Target ID" : "tophit", "#Query ID" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    df["tf_name"] = tf_name 

    # Grouping on smallest p-val or best p-pval for each motif:
    df["best_pval"] = df.groupby(["tf_name", "motif_id"])["E-value"].transform(min)
    least_pval_idx = df.groupby(["tf_name", "motif_id"])["E-value"].idxmin().values
    final_df = df.loc[least_pval_idx]
    select_cols = ["tophit", "motif_id", "Orientation", "p-value", "E-value"]
    tomtom_select_df = final_df.loc[:,select_cols]
    tomtom_select_df = tomtom_select_df.sort_values(["motif_id"]).reset_index(drop=True)
    
    return(tomtom_select_df)


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

    motif_ideas_final_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.bed"), header=True, index=True, sep="\t")
    motif_ideas_final_df.to_pickle(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.pkl"))
    motif_ideas_distribution.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_motifsfraction_ideas_piechart_dist.txt"), header=True, index=True, sep="\t")

    return(motif_ideas_final_df)


def generate_motifs_binned_coords(motifs_coordinates_info, tf_name, upstream_range, downstream_range, bin_size):
    #upstream_range = 2000
    #downstream_range = 2000
    #bin_size=50
    #motifs_coordinates_info = motifs_ideas_coord_df

    motifs_df =  motifs_coordinates_info.sort_values(["chrom","start","end"])
    upstream = upstream_range
    downstream = downstream_range
    bin_size = bin_size
    nrows =  motifs_df.shape[0]

    bins = range(-upstream, (downstream), bin_size)
    bin_len = len(bins)
    motifs_concat_df = pd.concat([motifs_df]*bin_len, ignore_index="TRUE")
    motifs_sorted_df = motifs_concat_df.sort_values(["chrom","start","end"])

    ### Copy the bin list that is deep copy:
    bin_start_list = bins[:]
    bin_end_list = []
    for each in bin_start_list:
        bin_end_list.append(each+bin_size)

    bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
    bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
    bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

    ### Combine the motifs df and bin df by cbind or column wise:
    temp_motifs_df = pd.concat([motifs_sorted_df.reset_index(), bin_concat_df], axis = 1)
    temp_motifs_df["motifs_midpoint"] = (temp_motifs_df["start"] + temp_motifs_df["end"])/2
    temp_motifs_df["motifs_midpoint"] = temp_motifs_df["motifs_midpoint"].round().astype(int)
    #final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "fimo_seq"]]
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "fimo_seq", "ideas_state"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]
    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'motifs_midpoint', u"motif_id", u"fimo_seq", "ideas_state"]
    final_motifs_df = final_motifs_df.loc[:,select_cols]

    ### Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    final_motifs_df = final_motifs_df.loc[final_motifs_df["chrom_start"] > 0, :] 
    final_motifs_df.to_csv(join(output_dir, tf_name + "_binned_motifs_coordinate_info_w_IDEAS.bed"), sep="\t", index=False, header=False)

    return(final_motifs_df)


############################
# Function implementations #
############################

fimo_master_tfdict = {}
motif_records_dict = {}
regex_pat = re.compile(r'.*fimo_motif_(.*).txt$')
# fimo_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL159319_SE_VS_SL159326_fimo_motif_FOXK1[FLAG].txt"
for fimo_file in fimo_filelist: #
    TF_name = regex_pat.findall(basename(fimo_file))[0] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
    bed_file = bed_filedict[TF_name]
    centrimo_file = centrimo_filedict[TF_name]
    tomtom_file = tomtom_filedict[TF_name]
    print("Currently processing : {} TF\n".format(TF_name))

    fimo_coord_df = final_fimo_motifs_model(fimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)       
    if TF_name not in fimo_master_tfdict:
        fimo_master_tfdict[TF_name] = fimo_coord_df
    
    # Increases search-space up-and-downstream of motifs:
    # fimo_coord_df["chromStart"] = fimo_coord_df["chromStart"] -1
    # fimo_coord_df["chromEnd"] = fimo_coord_df["chromEnd"] +1

    # Intersection of motif and total peak file:
    motif_bedfile = pybedtools.BedTool.from_dataframe(fimo_coord_df)
    peaks_df = pd.read_csv(bed_filedict[TF_name], sep="\t", header=None)
    peaks_df.rename(columns={0: "chr", 1: "start", 2: "end"}, inplace=True)
    sorted_peak_df = peaks_df.sort_values(["chr", "start", "end"])
    peaks_bedfile = pybedtools.BedTool.from_dataframe(sorted_peak_df)
    motif_peak_intersect = motif_bedfile.intersect(peaks_bedfile, wa=True, wb=True)
    
    # Read intersect file for further processing:
    motif_peak_df = pd.read_csv(motif_peak_intersect.fn, sep="\t", header=None)
    motif_peak_df.columns = fimo_coord_df.columns.tolist() + peaks_df.columns.tolist()

    motif_fraction_dict = {}        
    grouped_df = motif_peak_df.groupby(["motif_id"])
    gp_keylist = grouped_df.groups.keys()
    for motif_id in sorted(gp_keylist):
        motif_df = grouped_df.get_group(motif_id)
        motif_peak_df.groupby(["motif_id"]).size()
        motif_final_df = motif_df.drop_duplicates(["chr","start","end"]).reset_index(drop=True)
        motif_fraction_dict[motif_id] = [motif_final_df.shape[0], peaks_df.shape[0]]

    motif_dict_df = pd.DataFrame(motif_fraction_dict.items())
    motif_dict_df.columns = [ "motif_id", "peak_count" ]
    motif_dict_df[["peak_count", "total_peakcount"]] = motif_dict_df["peak_count"].apply(pd.Series) 
    motif_dict_df["fraction_bind"] = motif_dict_df["peak_count"]/motif_dict_df["total_peakcount"].astype(float)
    
    # Processed fimo files, ready for merge:
    fimo_df = motif_dict_df.sort_values(["motif_id"]).reset_index(drop=True)
    name = "fimo_"
    fimo_df = fimo_df.add_prefix(name)
    fimo_df.rename(columns={name+"motif_id": "motif_id"}, inplace=True) # for common-val based merge

    # Processing tomtom files for merge with fimo records:
    tomtom_df = final_tomtom_motifs_model(tomtom_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)
    name = "tom_"
    tomtom_df = tomtom_df.add_prefix(name)
    tomtom_df.rename(columns={name+"motif_id": "motif_id"}, inplace=True) # for common-val based merge

    # Processing centrimo files for merge with fimo records:
    centrimo_df = final_centrimo_motifs_model(centrimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)
    name = "cmo_"
    centrimo_df = centrimo_df.add_prefix(name)
    centrimo_df.rename(columns={name+"motif_id": "motif_id"}, inplace=True) # for common-val based merge

    # Merge fimo peak fraction, tomtom df and centrimo for each TF:
    dfs = [fimo_df, tomtom_df, centrimo_df]
    merged_ftc_df = reduce(lambda left,right: pd.merge(left,right, on='motif_id', how="outer"), dfs)
    # merged_df = pd.merge(fimo_tomtom_df, centrimo_df, on="motif_id", how="outer")

    # Insert final df to dictionary:
    motif_records_dict[TF_name] = merged_ftc_df

# Concat all TF's data to compile merged fimo-tomtom-centrimo df:       
merged_final_df = pd.concat(motif_records_dict).reset_index().drop(["level_1"],axis=1)
merged_final_df.rename(columns={"level_0":"TF_NAME", "motif_id": "MOTIF_ID"},inplace=True)
merged_final_df.to_csv(join(output_dir, "fimo_cisbp_tomtom_centrimo_df.txt"), sep="\t", header=True, index=False)
merged_final_df.to_excel(join(output_dir, "fimo_cisbp_tomtom_centrimo_df.xls"), header=True, sheet_name="Sheet1")
merged_final_df.to_pickle(join(output_dir,"fimo_cisbp_tomtom_centrimo_df.pkl"))
merged_final_df = pd.read_pickle(join(output_dir,"fimo_cisbp_tomtom_centrimo_df.pkl"))

# Merge meme motif sig and merged fimo-tomtom-centrimo:
combined_motifinfo_df = pd.merge(df_motif_table, 
                                merged_final_df, 
                                left_on=["TF_NAME", "MOTIF_ID"], right_on=["TF_NAME", "MOTIF_ID"], 
                                how="outer", 
                                indicator=True )
combined_motifinfo_df.to_csv(join(output_dir, "final_allcombined_motifinfo_cisbp_df.txt"), sep="\t", header=True, index=False)
combined_motifinfo_df.to_excel(join(output_dir, "final_allcombined_motifinfo_cisbp_df.xls"), header=True, sheet_name="Sheet1")
combined_motifinfo_df.to_pickle(join(output_dir,"final_allcombined_motifinfo_cisbp_df.pkl"))

# Filter significant motif info from all:
combined_motifinfo_df = pd.read_pickle(join(output_dir,"final_allcombined_motifinfo_cisbp_df.pkl"))
combined_significant_motifinfo_df = combined_motifinfo_df.loc[(combined_motifinfo_df["E-value"].astype(float) < 1e-05) \
                                    & (combined_motifinfo_df["cmo_E-value"].astype(float) < 1e-10) \
                                    & (combined_motifinfo_df["cmo_bin_width"].astype(float) < 150) \
                                    & (combined_motifinfo_df["SITES"].astype(float) >= 50)] # at least 10% of sites

combined_significant_motifinfo_df.drop(["SITES", "E-value_sig", "_merge", "cmo_altname","tom_p-value"], axis=1, inplace=True)
combined_significant_motifinfo_df.rename(columns={"fimo_total_peakcount" : "tf_peakcount"}, inplace=True)
combined_significant_motifinfo_df.rename(columns={"fimo_fraction_bind" : "motif_peak_fraction"}, inplace=True)
combined_significant_motifinfo_df["MOTIF_ID"].value_counts() # distribution of motifs
combined_significant_motifinfo_df.to_pickle(join(output_dir,"combined_significant_motifinfo_cisbp_df.pkl"))

# Filter motifs with DNA binding domains:
combined_significant_motifinfo_df = pd.read_pickle(join(output_dir,"combined_significant_motifinfo_cisbp_df.pkl"))
final_crcf_motifinfo = combined_significant_motifinfo_df.loc[~combined_significant_motifinfo_df["TF_NAME"].isin(dbf_anno_df["Target"])]
final_dbf_motifinfo = combined_significant_motifinfo_df.loc[combined_significant_motifinfo_df["TF_NAME"].isin(dbf_anno_df["Target"])]
final_dbf_motifinfo.to_csv(join(output_dir, "Final_dbf_motifinfo_cisbp.txt"), sep="\t", header=True, index=False)
final_dbf_motifinfo.to_excel(join(output_dir, "Final_dbf_motifinfo_cisbp.xls"), header=True, sheet_name="Sheet1")
final_dbf_motifinfo[final_dbf_motifinfo["tom_E-value"] > 0.05].to_excel(join(output_dir, "Final_dbf_motifinfo_cisbp_greaterthan_0.05similarity.xls"), header=True, sheet_name="Sheet1")

# For heatmap data:
final_dbf_motifinfo["annotation"] = final_dbf_motifinfo["TF_NAME"] + "." + final_dbf_motifinfo["MOTIF_ID"] 
final_dbf_motifinfo["tom_E-value"].fillna(1, inplace=True) # replace NA with 1 i.e no comparable motifs in DB
final_dbf_motifinfo_sorted = final_dbf_motifinfo.sort_values(["tom_E-value"]).reset_index(drop=True)

# You can call `quantile(i)` to get the i'th quantile,
# where `i` should be a fractional number.
final_dbf_motifinfo["tom_E-value"].quantile(0.1) # 10th percentile
final_dbf_motifinfo["tom_E-value"].quantile(0.5) # same as median
final_dbf_motifinfo["tom_E-value"].quantile(0.9) # 90th percentile

# For list of quantiles:
# Get quantiles for [.1, .2, .3, .4, .5, .6, .7, .8, .9]
import pylab as pl
percentile_arr = pl.frange(0.005,1,0.005)
tomtom_fraction_table_desc = final_dbf_motifinfo["tom_E-value"].describe(percentiles=percentile_arr)
tomtom_fraction_table = final_dbf_motifinfo["tom_E-value"].quantile(percentile_arr)
tomtom_fraction_df = pd.DataFrame(tomtom_fraction_table).reset_index()
tomtom_fraction_df.columns = ["motif_fraction", "tom_E-value"]
tomtom_fraction_df.to_csv(join(output_dir, "Novel_motifs_cisbp_tomtom_fraction_table_redo.txt"), sep="\t", header=True, index=False)
final_dbf_motifinfo["tom_E-value"].median() # 50 percentile
# final_dbf_motifinfo["tom_E-value"].mean()


##############################################################
# Pygrep to search the identifiers of best scoring motif hit:
##############################################################
#!/usr/bin/env python
 
# USAGE: 
# ./pygrep.py M4629_1.02 /gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme
 
# Using Loops for multiple strings:
# for each in M1925_1.02 M4629_1.02 M5323_1.02 M4481_1.02;do ./pygrep.py $each /gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme; done
 
import os,re, sys

file = open(sys.argv[2], "r")

tf_search_list = []
#for line in open(file).readlines():
for line in file:
    if re.search(sys.argv[1], line):
        print line,
        print line.split(),
        tf_search_list.append(line.split()[2])

print(tf_search_list)


#############################################

# For motif similarity scatter plot analyis:

#############################################


def final_tomtom_motifs_model_for_topScan(tomtom_file, tf_name, num_top_hit_scan, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)
    df = pd.read_csv(tomtom_file, sep="\t")
    df.rename(columns={"Target ID" : "tophit", "#Query ID" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    df["tf_name"] = tf_name 
    
    # Grouping on smallest p-val or best p-pval for each motif:
    # final_tomtom_df = df.set_index("tophit").groupby(["tf_name", "motif_id"])["E-value"].nsmallest(5).reset_index()
    df_list = []
    item_num = num_top_hit_scan
    select_cols = ["tf_name", "motif_id", "tophit", "Orientation", "E-value"]
    df = df.loc[:,select_cols]
    grp_object = df.groupby(["tf_name", "motif_id"])
    for key, val in grp_object:
        grp_df = grp_object.get_group(key)
        grp_df_sorted = grp_df.sort_values(["E-value"]).head(item_num)
        df_list.append(grp_df_sorted)
    combined_grp_df = pd.concat(df_list, ignore_index=True)
    
    return(combined_grp_df)

# Processing tomtom files for merge with fimo records:
tomtom_df_list = []
for key, val in tomtom_filedict.iteritems():
    print("\nProcessing {}".format(key))
    TF_name = key
    tomtom_file = val
    tomtom_df = final_tomtom_motifs_model_for_topScan(tomtom_file, TF_name, 50, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)
    name = "tom_"
    tomtom_df = tomtom_df.add_prefix(name)
    tomtom_df.rename(columns={name+"motif_id": "motif_id"}, inplace=True) # for common-val based merge
    tomtom_df_list.append(tomtom_df)

combined_tomtom_df = pd.concat(tomtom_df_list)

# Find identifier(tf_name) for each CisBP id:
#infile = open("/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens.meme", "rb")
infile = open("/gpfs/gpfs1/home/schhetri/Tools/final_custom_CisBP_database/Homosapiens_custom_cisbp_allmotifs.meme", "rb")
motif_list = []
for line in infile.readlines():
    if re.search("MOTIF", line):
        motif_list.append(line.rstrip())
        print(line)

cisbp_df = pd.Series(motif_list).str.split(expand=True)
cisbp_df.columns = ["motif", "cisbp_motif_id", "tf_name"]

# Annotate cisbp id from tophit tomtom df by merge with cisbp df:
tomtom_cisbp_df = pd.merge(combined_tomtom_df,cisbp_df, left_on="tom_tophit", right_on="cisbp_motif_id")
tomtom_cisbp_df.drop(["motif"], axis=1, inplace=True)
tomtom_cisbp_df_srt = tomtom_cisbp_df.sort_values(["tom_tf_name", "motif_id"]).reset_index(drop=True)

# Previous DBF 293 motifs:
#output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
combined_significant_motifinfo_df = pd.read_pickle(join(output_dir,"combined_significant_motifinfo_cisbp_df.pkl"))
final_crcf_motifinfo = combined_significant_motifinfo_df.loc[~combined_significant_motifinfo_df["TF_NAME"].isin(dbf_anno_df["Target"])]
final_dbf_motifinfo = combined_significant_motifinfo_df.loc[combined_significant_motifinfo_df["TF_NAME"].isin(dbf_anno_df["Target"])]

# Re-Annotate cisbp id from tophit tomtom df by merge with cisbp df:
final_dbf_cisbp_df = pd.merge(final_dbf_motifinfo, tomtom_cisbp_df_srt, left_on=["TF_NAME","MOTIF_ID"], right_on=["tom_tf_name", "motif_id"], how="left")
final_dbf_cisbp_df["TF_NAME_2"] = final_dbf_cisbp_df["TF_NAME"].str.replace("_human", "").str.replace("\[FLAG\]", "").str.upper().str.replace(r"_.*", "")
final_dbf_cisbp_df["tf_name_2"] = final_dbf_cisbp_df["tf_name"].str.upper().str.replace(r"_.*", "").str.replace("\(", "").str.replace("\)", "").str.replace(r"\..*", "")

true_df_list = []
false_df_list = []

dbf_grp_object = final_dbf_cisbp_df.groupby(["TF_NAME", "MOTIF_ID"])
for key, val in dbf_grp_object:
    print key, val
    grp_df = dbf_grp_object.get_group(key)
    if any(grp_df["TF_NAME_2"] == grp_df["tf_name_2"]):
        true_df = grp_df.loc[grp_df["TF_NAME_2"] == grp_df["tf_name_2"]]
        true_df_list.append(true_df.sort_values(["tom_E-value_y"]).head(1))
    else:
        false_df = grp_df.sort_values(["tom_E-value_y"]).head(1)
        false_df_list.append(false_df)


combined_true_df = pd.concat(true_df_list)
combined_true_df["anno"] = "match"
combined_false_df = pd.concat(false_df_list)
combined_false_df["anno"] = "mismatch"
combined_false_df.loc[combined_false_df["tom_E-value_y"].isnull(), "anno"] = "no_match"

combined_true_false_df = pd.concat([combined_true_df, combined_false_df], ignore_index=True)
combined_true_false_df_srt = combined_true_false_df.sort_values(["anno"])

# Order by tomtom E value for ordering within the group:
combined_grp_object = combined_true_false_df_srt.groupby(["anno"])
final_df_list = []
for key, val in combined_grp_object:
    final_df = combined_grp_object.get_group(key).sort_values(["tom_E-value_y"])
    final_df_list.append(final_df)

combined_tomtom_cisbp_scan_df = pd.concat(final_df_list,ignore_index=True)
combined_tomtom_cisbp_scan_df["MOTIF_ID_2"] = combined_tomtom_cisbp_scan_df["MOTIF_ID"].str.replace("MOTIF", "")
combined_tomtom_cisbp_scan_df["tf_motif_id"] = combined_tomtom_cisbp_scan_df["TF_NAME_2"] + "." + combined_tomtom_cisbp_scan_df["MOTIF_ID_2"]
combined_tomtom_cisbp_scan_df["tom_E-value_y"].fillna(1, inplace=True)
combined_tomtom_cisbp_scan_df["anno"] = combined_tomtom_cisbp_scan_df["anno"].replace({"match" : "concordant", "mismatch" : "discordant"})
combined_tomtom_cisbp_scan_df["tom_tophit_y"] = combined_tomtom_cisbp_scan_df["tom_tophit_y"].str.replace("@.*", "")
combined_tomtom_cisbp_scan_df.to_csv(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_redo_final.txt"), sep="\t", header=True, index=False)
combined_tomtom_cisbp_scan_df.to_excel(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_redo.xls"), header=True, index=False)

combined_tomtom_cisbp_scan_df.to_pickle(join(output_dir,"combined_tomtom_cisbp_scan_df_top50_redo_final.pkl"))
combined_tomtom_cisbp_scan_df = pd.read_pickle(join(output_dir, "combined_tomtom_cisbp_scan_df_top50_redo_final.pkl"))

combined_tomtom_cisbp_scan_df["anno"].value_counts()
combined_tomtom_cisbp_scan_df[combined_tomtom_cisbp_scan_df["anno"] == "concordant"]["TF_NAME_2"].nunique()
combined_tomtom_cisbp_scan_df[combined_tomtom_cisbp_scan_df["anno"] == "discordant"]["TF_NAME_2"].nunique()
combined_tomtom_cisbp_scan_df[combined_tomtom_cisbp_scan_df["anno"] == "no_match"]["TF_NAME_2"].nunique()
mismatch_df = combined_tomtom_cisbp_scan_df[combined_tomtom_cisbp_scan_df["anno"] == "mismatch"]
mismatch_df.sort_values(["tf_name_2"])

