import os, errno
import re
from glob import glob
from os.path import splitext, join, basename, dirname, expanduser
import pandas as pd, numpy as np
import pybedtools

# Overlap b/w Ryne's peak and our motif-based hotsites:
""" Local Machine analysis """
in_dir="/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"

df_ryne = pd.read_csv(join(in_dir, "MergedBedswCounts.txt"), sep="\t")
select_cols = df_ryne.columns.tolist()
[ select_cols.remove(each) for each in ["total_peaks", "TFs_Bound"]]
df_ryne = df_ryne.loc[:,select_cols]
header1 = df_ryne.columns.tolist() 

df_motif = pd.read_csv(join(in_dir, "Hotmotif_sites.bed"), sep="\t")
select_cols = [u'chrom', u'chromStart', u'chromEnd', u'uniq_tfcount']
df_motif = df_motif.loc[:,select_cols]
header2 = df_motif.columns.tolist() 

motif_pybed = pybedtools.BedTool.from_dataframe(df_motif)
ryne_pybed = pybedtools.BedTool.from_dataframe(df_ryne)
pybed_intersect = motif_pybed.intersect(ryne_pybed, wao=True)
df_intersect = pd.read_csv(pybed_intersect.fn, sep="\t", header=None)  
final_header = header2 + header1 + ["bp_overlap"]
df_intersect.columns = final_header
df_intersect_final = df_intersect.drop_duplicates([u'chrom', u'chromStart', u'chromEnd', u'uniq_tfcount']) 
df_intersect_final.equals(df_intersect)
df_hotsite_50 = df_intersect_final.loc[df_intersect_final["uniq_tfcount"]>=50]
df_hotsite_50_select = df_hotsite_50.loc[~(df_hotsite_50["bp_overlap"] == 0)]
df_hotsite_50_select["count_diff"] = df_hotsite_50_select["unique_peaks"].astype(int) - df_hotsite_50_select["uniq_tfcount"].astype(int)
df_hotsite_50_select.to_csv(join(in_dir, "Hotsites_motif_peak_based_50TF_more.bed"), sep="\t", header=True, index=False)
df_hotsite_50.to_csv(join(in_dir, "Hotsites_motif_based_50TF_more.bed"), sep="\t", header=True, index=False)

# For Chris - 100 or more TF bounds:
df_hotsite[df_hotsite_50["uniq_tfcount"] >=100]
df_hotsite_forChris = df_hotsite_50[df_hotsite_50["uniq_tfcount"] >=100].iloc[:,0:4]
df_hotsite_forChris.to_csv(join(in_dir, "HotMotifSites_100TFs_or_more_forChris.bed"), sep="\t", header=True, index=False)

# Confirming the overlap b/w df_hotsite_50 and df_overlap based on 50TF bounds:
df_ryne_50 = df_ryne.loc[df_ryne["unique_peaks"] >= 50]
df_motif_50 = df_motif.loc[df_motif["uniq_tfcount"] >= 50]
ryne_50_pybed = pybedtools.BedTool.from_dataframe(df_ryne_50)
motif_50_pybed = pybedtools.BedTool.from_dataframe(df_motif_50)
df_overlap = motif_50_pybed.intersect(ryne_50_pybed, u=True)


