#!/usr/bin/env python2.7

import os, errno, time
import re
from glob import glob
from os.path import splitext, join, basename, dirname, expanduser
import pandas as pd, numpy as np
import pybedtools
start_time = time.time()

# tf_files = expanduser("~/Dropbox/for_genemodels/idr_passed_peaks_total/test_analysis/SL*narrowPeak*")
tf_files = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*"

# genecode_file = expanduser("~/Dropbox/for_genemodels/gencode.v19.annotation.gtf")
genecode_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/for_genecodeTssProfile/gencode.v19.annotation.gtf"

# processing genecode(GTF) file:
genecode_file_edited = join(dirname(genecode_file), "gencode.v19.annotation_edit1.gtf")

if not os.path.exists(genecode_file_edited):
    os.system( ''' tail -n+6 %s | awk 'BEGIN{print "chrom\tstart\tend\tgene_name\tgene_id\tstrand\tfeature"} \
            {if($3=="gene")print $1,$4,$5,$18,$10,$7,$3}' OFS="\t" > %s ''' %(genecode_file, genecode_file_edited))

# output_dir = os.path.expanduser("~/Dropbox/for_genemodels")
output_dir = os.path.expanduser("~/for_chris/batch_I/for_genecodeTssProfile/for_genemodels")

suboutput_dir = join(output_dir, "tss_tf_intersect_files")

# download HepG2 FPKM file from encode on fly:
download_dir = output_dir
if not os.path.exists(join(download_dir,"ENCFF139ZPW.tsv")):
    os.system("wget -P {} https://www.encodeproject.org/files/ENCFF139ZPW/@@download/ENCFF139ZPW.tsv".format(download_dir))

# if not os.path.exists(output_dir):
#   os.makedirs(output_dir)

# Useful for job submission - if multiple threads running with a race condition to create the dir:
try:
    os.makedirs(suboutput_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


def generate_singlebinned_tss_coords_genemodel(genecode_edited_file, upstream_range, downstream_range):
    """ 
    Function to filter raw genecode model for strand specific coordinate settings 
    
    Args:
        data_raw(TSV bedfile): Raw genecodeV19 gtf formatted bed file with selective cols
    
    Returns:
        data(pd.DataFame): Processed and cleansed genecode model data and defined bin of TSS
    """

    # Read Tab-separated value bedfiles:
    df = pd.read_csv(genecode_file_edited, sep="\t")
    print("Current dimension of the genecode model:{}".format(df.shape))
    df = df.drop_duplicates()
    print("Dropping duplicates,if any - current dimension of genecode model:{}\n\n".format(df.shape))
    regex_pat = re.compile(r";", flags=0) # flags=re.IGNORECASE
    df["gene_name"] = df["gene_name"].str.replace(regex_pat, "")
    df["gene_id"] = df["gene_id"].str.replace(regex_pat, "")
    df_pos = df.loc[df["strand"] == "+"]; df_pos.head()
    df_neg = df.loc[df["strand"] == "-"]; df_neg.head()

    # Set up-and-downstream coordinates based on pos and neg strand
    df_pos["tss_midpoint"] = df_pos["start"]
    df_pos["chrom_Start"] = df_pos["tss_midpoint"].astype(int) - upstream_range # -3000 
    df_pos["chrom_End"] = df_pos["tss_midpoint"].astype(int) + downstream_range # +3000
    df_neg["tss_midpoint"] = df_neg["end"]
    df_neg["chrom_start"] = df_neg["tss_midpoint"].astype(int) + upstream_range # +3000
    df_neg["chrom_end"] = df_neg["tss_midpoint"].astype(int) - downstream_range # -3000
    
    # Maintain chrom_start coord > chrom_neg to avoid bed-intersection upset 
    df_neg.rename(columns={"chrom_start" : "chrom_End", "chrom_end" : "chrom_Start"}, inplace=True)
    
    # Combine positive and negative stranded genes:
    select_cols = [u'chrom', u'chrom_Start', u'chrom_End', u'start', u'end', u'strand', u'gene_name', u'gene_id', u'tss_midpoint']
    genecode_pos_df = df_pos.loc[:, select_cols]
    genecode_neg_df = df_neg.loc[:, select_cols]
    genecode_model_df = pd.concat([genecode_pos_df, genecode_neg_df]).sort_values(["chrom", "chrom_Start","chrom_End"])

    # Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    genecode_model_df.loc[genecode_model_df["chrom_Start"] < 0, "chrom_Start"] = 0
    tss_coord_df = genecode_model_df.copy()
    tss_coord_df.to_csv(join(output_dir, "genecodeV19_singlebinnned_tss_coordinate_info.bed"), sep="\t", index=False, header=False)
    
    return(tss_coord_df)


def load_combined_tf_pybedtool_object(tf_bedfile_df):
    """
        Function to process and generate pybedtool object for TFs bed file

        Returns:
            Pybedtool object with sorted bedfile
    """
    #tf_bedfile_df = pd.read_csv(tf_file_fullpath, sep="\t", header=None)
    tf_bedfile_df.rename(columns={0:"chrom", 1:"chromstart", 2:"chromend"}, inplace=True)
    ncols = tf_bedfile_df.shape[1]
    
    # tf_bedfile_df = tf_bedfile_df.iloc[:, [0,1,2]]
    # tf_bedfile_df.columns = ["chrom", "chromstart", "chromend"]
    sorted_df = tf_bedfile_df.sort_values(["chrom", "chromstart", "chromend"]).reset_index(drop=True) 
    sorted_df["peaks_midpoint"] = (sorted_df["chromstart"] + sorted_df["chromend"])/2
    sorted_df["peaks_midpoint"] = sorted_df["peaks_midpoint"].round().astype(int)
    sorted_df["start"] = sorted_df["peaks_midpoint"] - 50
    sorted_df["end"] = (sorted_df["peaks_midpoint"] + 50)
    sorted_df = sorted_df.drop_duplicates(["chrom", "start", "end", "tf_name"]) # if any
    select_cols = ["chrom","start", "end"] + range(3,ncols)
    sorted_df = sorted_df.loc[:,select_cols]
    tf_bedfile = pybedtools.BedTool.from_dataframe(sorted_df)

    return(tf_bedfile)


def generate_peaks_singlebinned_perc_bind(tssbins_coord_df, tfs_file_list, **kwargs):
    """ Function to generate fraction of peaks binding within defined region of Tss
        
        Args:
            tssbins_coord_df(pd.DataFrame): Genecode Tss binned dataframe
            tfs_file_list(pythonlist): List of TF peak files with full path
        
        Returns:
            Fraction of TF binding peaks within requested defined Tss region 
    """     
    print("kwargs: {}\n".format(kwargs)) # file_name =  kwargs["files_basename"]
    tss_sorted_df = tssbins_coord_df.sort_values(["chrom", "chrom_Start", "chrom_End"]).reset_index(drop=True)
    tssbins_bedfile = pybedtools.BedTool.from_dataframe(tss_sorted_df)

    master_dict = {}
    for idx, each_file in enumerate(tfs_file_list):
        regex_pat = re.compile(r".*narrowPeak_")
        file_basename = regex_pat.split(basename(each_file))[-1]

        # Calling other function to generate pybedtool object and do intersection:
        tf_bedfile = load_tf_pybedtool_object(each_file)
        print("Processing the intersection for: {}\n".format(file_basename))
        pybed_outfile = join(suboutput_dir, (file_basename + "_genecodeTss_intersect.bed"))
        pybed_outfile_v = join(suboutput_dir, (file_basename + "_genecodeTss_outersect.bed"))

        # if not os.path.exists(pybed_outfile):
        tss_tf_intersect = tssbins_bedfile.intersect(tf_bedfile, wa = True, wb = True, output=pybed_outfile)        
        print(tss_tf_intersect.head())  
        # tf_bedfile.intersect(peaks_bedfile, wa = True, wb = True, v = True, output=pybed_outfile_v)

        # Working with the dataframes; reading output file of pybedtool intersect:
        tss_tf_df = pd.read_csv(tss_tf_intersect.fn, sep = "\t", header = None)
        final_df = tss_tf_df.iloc[:,[6,7,8]].drop_duplicates().reset_index(drop=True)
        
        # For unmerged TSS dataset, uncomment below:
        #final_df = tss_tf_df.iloc[:,[8,9,10]].drop_duplicates().reset_index(drop=True)
        final_df.columns = ["chrom", "start", "end"]
        bind_count = final_df.shape[0]
        total_count = tf_bedfile.count()
        fraction_bind = round(bind_count/float(total_count), 4)
        print("Dimension of currently intersected peak file is: {}".format(final_df.shape))

        if file_basename not in master_dict:
          master_dict[file_basename] = [fraction_bind, total_count]
        print("Intersection of {} completed!!!...\n\n".format(file_basename))
    
    # Combine all dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    fraction_bind_combined_df = pd.DataFrame(master_dict.items())
    fraction_bind_combined_df.columns = ["tf_name", "bind_fraction"]
    fraction_bind_combined_df[["bind_fraction", "peak_count"]] = fraction_bind_combined_df["bind_fraction"].apply(pd.Series)
    fraction_bind_combined_df["peak_count"] = fraction_bind_combined_df["peak_count"].astype(int)
    print("Combined dataframe shape : {}".format(fraction_bind_combined_df.shape))
    fraction_bind_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = False)
    
    return(fraction_bind_combined_df)


def final_gene_expression_model(gene_expression_file):
    """ Function to filter raw gene expression file for downstream analysis 
    Args:
        data_raw(TSV bedfile): Raw gene_expression ENCODE bed file with selective cols
    Returns:
        data(pd.DataFame): Processed and cleansed gene expression model data
    """
    expression_df = pd.read_csv(gene_expression_file, sep="\t")
    expression_df["geneID_ensemble"] = expression_df["gene_id"].str.replace(r"\..*", ""); expression_df.head()
    return(expression_df)


def concat_tf_files(tfs_file_list):
    # Concat TFs file # tfs_filelist = glob(tf_files)
    concat_list = []
    for each_file in tfs_filelist:
        tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
        print("processing : {}\n".format(tf_name))
        tf_df = pd.read_csv(each_file, sep="\t", header=None)
        tf_df["tf_name"] = tf_name
        concat_list.append(tf_df)

    combined_tf_df = pd.concat(concat_list, ignore_index=True)
    combined_tf_df.to_csv(join(output_dir,"combined_tf_df.bed"), sep="\t", header=None, index=False)
    return(combined_tf_df)


# For fraction bind with expressed genes:
def main():
    """
        Function main to call sub-functions and automate all
    """
    tfs_filelist = glob(tf_files) # file list from TF filepath regex
    tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 2000, 2000)
    
    # Gene expression data with 1TPM or more:
    gene_exp_df = final_gene_expression_model(join(download_dir,"ENCFF139ZPW.tsv"))
    gene_exp_df = gene_exp_df.loc[:,["gene_id", "TPM", "FPKM"]]
    merged_df = pd.merge(tss_coord_df, gene_exp_df, on="gene_id")
    tss_coord_df = merged_df.loc[merged_df["TPM"]>=1]

    tss_pybed = pybedtools.BedTool.from_dataframe(tss_coord_df)
    merged_tss_pybed = tss_pybed.merge(c=[6,7,7,8,9,10,11], o=["collapse","count_distinct","collapse", \
                        "collapse","collapse","mean","mean"])
    merged_tss_coord_df = pd.read_csv(merged_tss_pybed.fn, sep="\t", header=None)
    merged_tss_coord_df.columns = ["chrom", "chrom_Start", "chrom_End", "strand", "gene_counts", "gene_name", "gene_id", "txStart", "TPM", "FPKM"]

    combined_tf_df = concat_tf_files(tfs_filelist)
    combined_tf_df = combined_tf_df.rename(columns={0:"chr", 1:"start", 2:"end"})
    
    # +/- 50bp from peak center:
    combined_tf_df["peaks_midpoint"] = (combined_tf_df["start"] + combined_tf_df["end"])/2
    combined_tf_df["peaks_midpoint"] = combined_tf_df["peaks_midpoint"].round().astype(int)
    combined_tf_df["start"] = combined_tf_df["peaks_midpoint"] - 50
    combined_tf_df["end"] = (combined_tf_df["peaks_midpoint"] + 50)

    combined_tf_df_sort = combined_tf_df.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    select_cols = ["chr", "start", "end", "tf_name"]
    sorted_tf_df = combined_tf_df_sort.loc[:,select_cols]
    sorted_tf_df["tf_name"] = sorted_tf_df["tf_name"].str.replace("\\[FLAG\\]", "").str.replace("_human", "")

    sorted_tf_bedfile = pybedtools.BedTool.from_dataframe(sorted_tf_df)
    #final_tss_bedfile = pybedtools.BedTool.from_dataframe(tss_coord_df)
    final_tss_bedfile = pybedtools.BedTool.from_dataframe(merged_tss_coord_df)
    tss_tf_intersect = final_tss_bedfile.intersect(sorted_tf_bedfile, wao = True)        
    tss_tf_df = pd.read_csv(tss_tf_intersect.fn, sep = "\t", header = None)
    tss_tf_df.columns = merged_tss_coord_df.columns.tolist() + sorted_tf_df.columns.tolist() + ["overlap"]
    print(tss_tf_df.head()) 

    # Tss with 0 or more than 0 peak overlap
    tss_tf_df_0 = tss_tf_df.loc[(tss_tf_df["overlap"] == 0)]; tss_tf_df_0.shape
    tss_tf_df_1 = tss_tf_df.loc[~(tss_tf_df["overlap"] == 0)]; tss_tf_df_1.shape

    # Groupby and generate count and uniq tfs:
    tss_tf_df_gpby = tss_tf_df_1.groupby(["chrom", "chrom_Start", "chrom_End", "txStart", "gene_name", "gene_id", "TPM", "FPKM"]) \
                    ["tf_name"].apply(lambda x: (x.unique(), x.nunique()))
    tss_tf_df_gpby_1 = tss_tf_df_gpby.reset_index(); tss_tf_df_gpby_1.shape
    tss_tf_df_gpby_1[["tf_name", "uniq_tf_count"]] = tss_tf_df_gpby_1["tf_name"].apply(pd.Series)

    # Merge 0 and more than 1 tf count containing dataframe:
    select_cols = tss_tf_df_gpby_1.columns.tolist()
    tss_tf_df_0 = tss_tf_df_0.loc[:,select_cols]
    tss_tf_df_0["uniq_tf_count"] = 0

    final_combined_df = pd.concat([tss_tf_df_gpby_1, tss_tf_df_0],  ignore_index=True)
    final_combined_df.shape
    final_combined_df.rename(columns={"tf_name": "tfs_associated"}, inplace=True)

    final_combined_df.sample(57).nlargest(5,"uniq_tf_count")

    return(final_combined_df)

if __name__ == '__main__':
    fraction_bind_data = main()



