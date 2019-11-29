import re, os 
from os.path import join, basename, expanduser
from glob import glob
import pandas as pd, numpy as np
import pybedtools 
from pyfasta import Fasta
from os import makedirs, rmdir, remove
from os.path import expanduser, exists

# Generate dict for peak-bed file list :
file_pat ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"
bed_filelist = glob(file_pat)
bed_regex_pat = re.compile(r'unique_TFs/.*narrowPeak_(.*)') 
bed_filedict = {bed_regex_pat.findall(file)[0]:file for file in bed_filelist}

# Generate dict for fimo file list :
file_pat = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/fimo_motifs_total/idr_passed_motifs_total/unique_TFs/SL*fimo_motif*"
fimo_filelist = glob(file_pat)
fimo_regex_pat = re.compile(r'.*unique_TFs/.*_fimo_motif_(.*).txt')
fimo_filedict = {fimo_regex_pat.findall(file)[0]:file for file in fimo_filelist}


def final_fimo_motifs_model(fimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["motif_mid_point"] = (df["chromStart"] + df["chromEnd"])/2 
	df["motif_mid_point"] = df["motif_mid_point"].astype(int)
	df["motif_halfsize"] = (df["chromEnd"] - df["chromStart"])/2
	df["motif_halfsize"] = df["motif_halfsize"].astype(int)
    df["tf_name"] = tf_name 
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_mid_point", "motif_halfsize", "motif_id", "tf_name", "strand", "p-value"]
    motif_select_df = df.loc[:,select_cols]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)

    return(motif_sorted_df)


def final_fimo_motifs_model_with_pval(fimo_file, tf_name, pval, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["chromStart"] = df["start"] + 1 # confirmed with pybed fasta to map right sequence 
    df["chromEnd"] =  df["stop"] + 2 # as bed intersection is exclusive of end coord
    df["motif_mid_point"] = (df["chromStart"] + df["chromEnd"])/2 
	df["motif_mid_point"] = df["motif_mid_point"].astype(int)
	df["motif_halfsize"] = (df["chromEnd"] - df["chromStart"])/2
	df["motif_halfsize"] = df["motif_halfsize"].astype(int)
    df["tf_name"] = tf_name 
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "chromStart", "chromEnd", "motif_mid_point", "motif_halfsize", "motif_id", "tf_name", "strand", "p-value"]
    motif_select_df = df.loc[:,select_cols]
    #motif_select_df = motif_select_df.loc[motif_select_df["p-value"] < pval]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "chromStart","chromEnd"]).reset_index(drop=True)

    return(motif_sorted_df)


def combine_tf_bedfiles_midpt_centered(TF_file, TF_name):
	tf_name = TF_name
	df_read = pd.read_csv(TF_file, header=None, sep="\t")
	df_read = df_read.iloc[:,[0,1,2]]
	df_read.columns = ["peak_chrom", "start", "end"]
	df_read["mid_point"] = (df_read["start"] + df_read["end"])/2 
	df_read["mid_point"] = df_read["mid_point"].astype(int)
	df_read["tf_name"] = tf_name
	select_cols = ["peak_chrom","peak_start","peak_end","start","end", "strand", "tf_name"]
	select_cols = ["peak_chrom","start","end", "mid_point", "tf_name"]
	df_read = df_read.loc[:,select_cols]
	df_final = df_read.sort_values(["peak_chrom","start","end"]).reset_index(drop=True)
	print df_final
	return(df_final)


############################
# Function implementations #
############################
# Can avoid the loop below using pickle load
# only file containing 293 motifs with top50 scan required:

output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis_tomtom_redo"
output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

empty_fimo_coord_dict = {}
fimo_master_tfdict = {}
for tf in fimo_filedict:
    TF_name = tf
    print("Currently processing : {} TF\n".format(TF_name))
    fimo_file = fimo_filedict[tf]
    #fimo_file = fimo_filedict["CTCF"]
    #TF_name = "CTCF"
    #tf = "CTCF"
    fimo_file = fimo_filedict[tf]
    #fimo_coord_df = final_fimo_motifs_model_with_pval(fimo_file, TF_name, 1e-05, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)       
    fimo_coord_df = final_fimo_motifs_model(fimo_file, TF_name, motif1=1, motif2=2, motif3=3, motif4=4, motif5=5)       

    if not fimo_coord_df.empty:
	    # Intersection of motif and total peak file:
	    motif_bedfile = pybedtools.BedTool.from_dataframe(fimo_coord_df)
	    peaks_df = combine_tf_bedfiles_midpt_centered(bed_filedict[TF_name], TF_name)
	    peaks_bedfile = pybedtools.BedTool.from_dataframe(peaks_df)
	    motif_peak_intersect = motif_bedfile.intersect(peaks_bedfile, wa=True, wb=True)
	    
	    # Read intersect file for further processing:
	    motif_peak_df = pd.read_csv(motif_peak_intersect.fn, sep="\t", header=None)
	    motif_peak_df.columns = fimo_coord_df.columns.tolist() + peaks_df.columns.tolist()
	    motif_peak_df["offset_dist"] = motif_peak_df["mid_point"] - motif_peak_df["motif_mid_point"]
		motif_peak_df["abs_offset_dist"] = abs(motif_peak_df["offset_dist"])

		# Find nearest motif to the peak
		gpd_df = motif_peak_df.groupby(["peak_chrom", "start", "end"])
		gpd_df_list = []
		for key,val in gpd_df:
			motif_peak_gp = gpd_df.get_group(key)
			motif_peak_gp = motif_peak_gp.nsmallest(1, "abs_offset_dist")
			gpd_df_list.append(motif_peak_gp)
			
		motif_peak_gp_df = pd.concat(gpd_df_list, ignore_index=True)
		motif_offset_df = motif_peak_gp_df.groupby(["motif_id"])["abs_offset_dist"].median().reset_index()
		motif_size_df = motif_peak_gp_df.groupby(["motif_id"])["motif_halfsize"].median().reset_index()
		merged_motif_df = pd.merge(motif_offset_df, motif_size_df, on=["motif_id"])
		merged_motif_df["adj_offset_dist"] = merged_motif_df["abs_offset_dist"] - merged_motif_df["motif_halfsize"]
		merged_motif_df.loc[merged_motif_df["adj_offset_dist"] < 0, "adj_offset_dist"] = 0

		if TF_name not in fimo_master_tfdict:
			fimo_master_tfdict[TF_name] = merged_motif_df
	else:
		print("\n{} has no significant fimo coordinates".format(TF_name))
		empty_fimo_coord_dict[TF_name] = np.nan
		continue

# Concat all TF's data to compile merged fimo-tomtom-centrimo df:       
combined_motif_df = pd.concat(fimo_master_tfdict).reset_index().drop(["level_1"],axis=1)
combined_motif_df.rename(columns={"level_0":"TF_NAME", "motif_id" : "MOTIF_ID"},inplace=True)
combined_motif_df["TF_NAME_2"] = combined_motif_df["TF_NAME"].str.replace("_human", "").str.replace("\[FLAG\]", "").str.upper().str.replace(r"_.*", "")
combined_motif_df["MOTIF_ID_2"] = combined_motif_df["MOTIF_ID"].str.replace("MOTIF", "")
combined_motif_df["tf_motif_id"] = combined_motif_df["TF_NAME_2"] + "." + combined_motif_df["MOTIF_ID_2"]
combined_motif_df.to_pickle(join(output_dir,"combined_motif_offset_peak_df_nearest.pkl"))
combined_motif_df = pd.read_pickle(join(output_dir,"combined_motif_offset_peak_df_nearest.pkl"))

# Import novel motif df: 
novel_filtered_motif_df = pd.read_csv(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_redo.txt"), sep="\t")
novel_filtered_motif_df.loc[novel_filtered_motif_df["TF_NAME"] == "MBD1_v1[FLAG]", "tf_motif_id"] = "MBD1_v1.1"
novel_filtered_motif_df.loc[novel_filtered_motif_df["TF_NAME"] == "MBD1_v2[FLAG]", "tf_motif_id"] = "MBD1_v2.1"

offset_motif_df = combined_motif_df[["tf_motif_id", "abs_offset_dist", "motif_halfsize", "adj_offset_dist"]]
final_merged_df = pd.merge(novel_filtered_motif_df, offset_motif_df, on=["tf_motif_id"])
final_merged_df["anno_new"] = final_merged_df["anno"].replace({"match":"concordance", "mismatch":"discordance"})
final_merged_df.to_csv(join(output_dir, "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset_nearest_redo.txt"), sep="\t", header=True, index=False)

print "Job completed successfully!!"





