import pandas as pd, numpy as np
import pybedtools
import pickle, scipy
import os, re
from glob import glob
from os.path import join, basename, splitext

main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/gatad2a_cobind_analysis"
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

""" All TF containing dir """
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob(dir_path) # represents all tf file list

gatad2a_peakpath = glob(dir_path + "GATAD2A*")[0]
gatad2a_cobind_list = ["ZNF219", "SMAD4", "RARA", "ARID5B", "FOXA3", "SOX13_iso1"]
gatad2a_cobind_peakspath = [ glob(dir_path + each + "*")[0] for each in gatad2a_cobind_list ]

def str_join(df, sep, *cols):
    from functools import reduce
    return reduce(lambda x, y: x.astype(str).str.cat(y.astype(str), sep=sep), [df[col] for col in cols])


def pybedfile_format(file_path, tf_name):
    peak_df = pd.read_csv(file_path, sep="\t", header=None)
    peak_df = peak_df.iloc[:, 0:3]
    peak_df.columns = ["chrom", "chromStart", "chromEnd"]
    peak_df["tf_name"] = tf_name
    peak_df = peak_df.sort_values(["chrom", "chromStart", "chromEnd"])
    bedfile = pybedtools.BedTool.from_dataframe(peak_df)
    return(bedfile)

gatad2a_bedfile = pybedfile_format(gatad2a_peakpath, "GATAD2A")

intersect_df_list = []
for each_tf in gatad2a_cobind_list:
    print("\nProcessing : {}".format(each_tf))
    cobind_tfpath = glob(dir_path + each_tf + "*")[0]
    cobind_bedfile = pybedfile_format(cobind_tfpath, each_tf)
    peaks_intersect = gatad2a_bedfile.intersect(cobind_bedfile, wao=True)
    # Read file and transform to binary format (0: nooverlap, 1: overlap)
    peaks_intersect_df = pd.read_csv(peaks_intersect.fn, sep="\t", header=None)
    peaks_intersect_df.columns = ["chrom", "chromStart", "chromEnd", "tf_name", \
                                    "chr", "start", "end", "tf_name2" ] + [each_tf]
    peaks_intersect_uniqdf = peaks_intersect_df.drop_duplicates(["chrom", "chromStart", "chromEnd"])
    peaks_intersect_uniqdf[each_tf] = np.where(peaks_intersect_uniqdf[each_tf] == 0, 0, 1)
    select_cols = ["chrom", "chromStart", "chromEnd"] + [each_tf]
    peaks_intersect_uniqdf = peaks_intersect_uniqdf.loc[:,select_cols]
    intersect_df_list.append(peaks_intersect_uniqdf)
    # peaks_intersect_uniqdf.loc[peaks_intersect_uniqdf[each_tf] > 0, each_tf ] = 1

merged_df = reduce(lambda left, right: pd.merge(left, right, on=["chrom","chromStart","chromEnd"]), intersect_df_list)
merged_df["Loci"] = str_join(merged_df, "_", "chrom", "chromStart", "chromEnd")
merged_df["GATAD2A"] = 1
select_cols = ["Loci", "GATAD2A"] + gatad2a_cobind_list
final_df = merged_df.loc[:,select_cols]
final_df.to_csv(join(main_dir, "merged_gatad2a_loci_for_coboundTF_overlap.txt"), sep="\t", header=True, index=False)

