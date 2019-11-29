import re, os 
from os.path import join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists

# import matplotlib
# import matplotlib.pyplot as plt
# from ggplot import *
# from plotnine import *
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
tomtom_filelist = glob(file_pat)
tomtom_regex_pat = re.compile(r'unique_TFs/.*_narrowPeak_(.*?)/meme_*')
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
    df["best_pval"] = df.groupby(["tf_name", "motif_id"])["p-value"].transform(min)
    least_pval_idx = df.groupby(["tf_name", "motif_id"])["p-value"].idxmin().values
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

    # motif_ideas_final_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.bed"), header=True, index=True, sep="\t")
    # motif_ideas_final_df.to_pickle(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.pkl"))
    # motif_ideas_distribution.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_motifsfraction_ideas_piechart_dist.txt"), header=True, index=True, sep="\t")

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


output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

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
merged_final_df.to_csv(join(output_dir, "fimo_tomtom_centrimo_df.txt"), sep="\t", header=True, index=False)
merged_final_df.to_excel(join(output_dir, "fimo_tomtom_centrimo_df.xls"), header=True, sheet_name="Sheet1")
merged_final_df.to_pickle(join(output_dir,"fimo_tomtom_centrimo_df.pkl"))
merged_final_df = pd.read_pickle("fimo_tomtom_centrimo_df.pkl")

# Merge meme motif sig and merged fimo-tomtom-centrimo:
combined_motifinfo_df = pd.merge(df_motif_table, 
                                merged_final_df, 
                                left_on=["TF_NAME", "MOTIF_ID"], right_on=["TF_NAME", "MOTIF_ID"], 
                                how="outer", 
                                indicator=True )
combined_motifinfo_df.to_csv(join(output_dir, "final_allcombined_motifinfo_df.txt"), sep="\t", header=True, index=False)
combined_motifinfo_df.to_excel(join(output_dir, "final_allcombined_motifinfo_df.xls"), header=True, sheet_name="Sheet1")
combined_motifinfo_df.to_pickle(join(output_dir,"final_allcombined_motifinfo_df.pkl"))
combined_motifinfo_df = pd.read_pickle(join(output_dir,"final_allcombined_motifinfo_df.pkl"))

# Concat all TF's data to compile fimo combined df:     
alltf_fimo_df = pd.concat(fimo_master_tfdict).reset_index().drop(["level_1"],axis=1)
alltf_fimo_df.rename(columns={"level_0":"TF_NAME", "motif_id": "MOTIF_ID"},inplace=True)
alltf_fimo_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_df.bed"), sep="\t", header=True, index=False)
alltf_fimo_df.to_pickle(join(output_dir,"final_tfcombined_fimocoords_df.pkl"))
alltf_fimo_df_merged = pd.merge(alltf_fimo_df, 
                    merged_final_df[["TF_NAME", "MOTIF_ID","fimo_fraction_bind", "fimo_total_peakcount"]],
                    left_on = ["TF_NAME", "MOTIF_ID"], right_on=["TF_NAME", "MOTIF_ID"],
                    how = "left") # unsorted df

# Sort fimo coords and drop extra TF_NAME col:
sorted_combinedfimo_df = alltf_fimo_df_merged.sort_values(["chrom", "chromStart", "chromEnd"]).reset_index(drop=True)
sorted_combinedfimo_df.drop(["TF_NAME"], axis=1, inplace=True) # getting rid of redundancy
sorted_combinedfimo_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind.bed"), sep="\t", header=True, index=False)
sorted_combinedfimo_df.to_pickle(join(output_dir,"final_tfcombined_fimocoords_with_fractionbind.pkl"))
sorted_combinedfimo_df = pd.read_pickle(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind.pkl"))

# Assign IDEAS annotation to sorted and all tf combined fimo coords:
retain_cols = [0,1,2,3,4,5,6,7]
ideas_combinedfimo_df = assign_IDEAS_State(sorted_combinedfimo_df, retain_cols)
ideas_combinedfimo_df.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_fractionbind_ideastates.bed"), sep="\t", header=True, index=False)
ideas_combinedfimo_df.to_pickle(join(output_dir,"final_tfcombined_fimocoords_with_fractionbind_ideastates.pkl"))

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

# ideas_head = ideas_combinedfimo_df.head()
ideas_combinedfimo_df["ideas_anno"] = ideas_combinedfimo_df["ideas_state"].replace(replace_ideas)
ideas_combinedfimo_df["ideas_reanno"] = ideas_combinedfimo_df["ideas_anno"].replace(remap_ideas)
ideas_combinedfimo_df["sp_motifcount"] = ideas_combinedfimo_df.groupby(["tf_name", "MOTIF_ID"]).transform(sum)

# Distribution/bindfraction across cis-regions - Groupby Tf name and motif id:
combinedfimo_dist = ideas_combinedfimo_df.groupby(["tf_name", "MOTIF_ID"]).size().reset_index(name="motifcount")
combinedfimo_dist_ideas = ideas_combinedfimo_df.groupby(["tf_name", "MOTIF_ID","ideas_anno"]).size().reset_index(name="motifcount_ideas")
combinedfimo_dist_ideas_merged = pd.merge(combinedfimo_dist, combinedfimo_dist_ideas, on=["tf_name", "MOTIF_ID"], how="left")
combinedfimo_dist_ideas_merged["cisbind_fraction"] = combinedfimo_dist_ideas_merged["motifcount_ideas"]/combinedfimo_dist_ideas_merged["motifcount"].astype(float)
combinedfimo_dist_ideas_merged.to_csv(join(output_dir, "final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t", header=True, index=False)
combinedfimo_dist_ideas_merged.to_pickle(join(output_dir,"final_tfcombined_fimocoords_with_ideascisreg_fractionbind.pkl"))

for key,values in combinedfimo_dist_ideas_merged.groupby(["MOTIF_ID"]):
    values.to_csv(join(output_dir, key + "_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t", header=True, index=False)
    heatmapdf = values.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
    heatmapdf.to_csv(join(output_dir, key + "_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=False)
    print key, values, heatmapdf

# Check Instance:
# foxk1_dist = combinedfimo_dist.loc[combinedfimo_dist["tf_name"] == "FOXK1[FLAG]"]
# foxk1_ideasdist = combinedfimo_dist_ideas.loc[combinedfimo_dist_ideas["tf_name"] == "FOXK1[FLAG]"]
# foxk1_merged = pd.merge(foxk1_ideasdist, foxk1_dist, on=["tf_name", "MOTIF_ID"], how="left")
# foxk1_merged["cisbind_fraction"] = foxk1_merged["motifcount_ideas"]/foxk1_merged["motifcount"].astype(float)
# foxk1_merged.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")

# Piechart based dist:
ideas_combinedfimo_dist = ideas_combinedfimo_df["ideas_anno"].value_counts().reset_index()
ideas_combinedfimo_dist.reset_index()["ideas_anno"].sum() == sorted_combinedfimo_df.shape[0]


##############################################
# Motif associated analysis and plots       ##
##############################################


# Fraction of Peaks accounted by each Motif:
# (Expected : In concordance with order of significantly associated motifs)

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis")
motif_df = pd.read_csv(join(output_dir, "final_allcombined_motifinfo_df.txt"), sep="\t")

# Plot 1:
motif_df_heat = motif_df.pivot(index="TF_NAME", columns="MOTIF_ID", values="fimo_fraction_bind")
motif_df_heat.to_csv(join(output_dir, "final_fimomotif_peakfraction_heatmapdata.txt"), sep="\t", header=True, index=True)

motif1_df = pd.read_csv(join(output_dir, "MOTIF1_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t")
motif1_df_heat = motif1_df.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
motif1_df_heat.to_csv(join(output_dir, "MOTIF1_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=True)


motif2_df = pd.read_csv(join(output_dir, "MOTIF2_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t")
motif2_df_heat = motif2_df.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
motif2_df_heat.to_csv(join(output_dir, "MOTIF2_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=True)


motif3_df = pd.read_csv(join(output_dir, "MOTIF3_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t")
motif3_df_heat = motif3_df.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
motif3_df_heat.to_csv(join(output_dir, "MOTIF3_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=True)


motif4_df = pd.read_csv(join(output_dir, "MOTIF4_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t")
motif4_df_heat = motif4_df.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
motif4_df_heat.to_csv(join(output_dir, "MOTIF4_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=True)


motif5_df = pd.read_csv(join(output_dir, "MOTIF5_final_tfcombined_fimocoords_with_ideascisreg_fractionbind.txt"), sep="\t")
motif5_df_heat = motif5_df.pivot(index="tf_name", columns="ideas_anno", values="cisbind_fraction")
motif5_df_heat.to_csv(join(output_dir, "MOTIF5_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=True, index=True)


#########################################################
# Common-Motifs available to MEME-cisBP-JASPAR centrimo #
#########################################################


test=$(find finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/*narrowPeak*/meme_cisbp_jaspar_centrimo_dir -type f -name "centrimo.txt")
for each in $test; do 
    base_dir=$(dirname $(dirname $each)); 
    find ${base_dir}/cisbp_singlefactor_based_meme -type f -name "*.meme"; 
done >> test.meme

cat $(cat test.meme) > concat.meme
grep "MOTIF" concat.meme

file_pat="finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis/*narrowPeak*/meme_cisbp_jaspar_centrimo_dir/centrimo.txt"
file_list=glob(file_pat) # 
tf_dir_list = []
for each in file_list:
    base_dir = dirname(dirname(each))
    tf_dir_list.append(base_dir)

cisbp_memefile_list = []
jaspar_memefile_list = []
for each_tf_dir in tf_dir_list:
    cisbp_memefile = glob(join(each_tf_dir,"cisbp_singlefactor_based_meme","*.meme"))
    jaspar_memefile = glob(join(each_tf_dir,"jaspar_singlefactor_based_meme","*.meme"))
    cisbp_memefile_list.append(cisbp_memefile)    
    jaspar_memefile_list.append(jaspar_memefile)    


cisbp_memefile = "./cisbp_motif_file_97factors.meme"


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

# Generate dict for peak bed file list :
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"
bed_filelist = glob(file_pat)
bed_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_filedict = {bed_regex_pat.findall(file)[0]:file for file in bed_filelist}


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

# fimo_master_tfdict = {}
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
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm/random5000_samples/gkm_train_output/*.cvpred.txt" # at motif level
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/kmer_svm_peaklevel/random5000_samples/gkm_train_output/*.cvpred.txt" # at full peak level
crossval_filelist = glob(file_pat)

regex_pat = re.compile(r'^(.*)_motifs_sample.cvpred.txt') # at motif level
regex_pat = re.compile(r'^(.*)_peaks_sample.cvpred.txt') # at full peak level
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
prauc_anno_df.to_csv(join(output_dir, "PRAUC_all_TF_annotated_peaklevel.txt"), sep="\t", header=True, index=False)

from plotnine import *
import pandas as pd
from os.path import join
""" Local machine plotting """

plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/PRAUC_all_TF_annotated.txt", sep="\t")
plot_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/PRAUC_all_TF_annotated_peaklevel.txt", sep="\t")

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

# # Same plot, but with R base plot:
# output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment"
# pdf(file.path(output_dir, "Mean_PRAUC_peaklevel_10X_nullsets.pdf"))
# plot_df <- fread("/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/PRAUC_all_TF_annotated_peaklevel.txt", sep="\t")
# boxplot(plot_df$pr_auc ~ plot_df$Category, main="Mean PR-AUC=0.66", col=c("blue", "green"))
# dev.off()

 
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

