

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


def load_tf_pybedtool_object(tf_file_fullpath):
    """
        Function to process and generate pybedtool object for TFs bed file

        Returns:
            Pybedtool object with sorted bedfile
    """
    regex_pat = re.compile(r".*narrowPeak_")
    file_basename = regex_pat.split(basename(tf_file_fullpath))[-1]
    print("Processing TF bed file: {}\n".format(file_basename))
    tf_bedfile_df = pd.read_csv(tf_file_fullpath, sep="\t", header=None)
    tf_bedfile_df.rename(columns={0:"chrom", 1:"chromstart", 2:"chromend"}, inplace=True)
    ncols = tf_bedfile_df.shape[1]
    
    # tf_bedfile_df = tf_bedfile_df.iloc[:, [0,1,2]]
    # tf_bedfile_df.columns = ["chrom", "chromstart", "chromend"]
    sorted_df = tf_bedfile_df.sort_values(["chrom", "chromstart", "chromend"]).reset_index(drop=True) 
    sorted_df["peaks_midpoint"] = (sorted_df["chromstart"] + sorted_df["chromend"])/2
    sorted_df["peaks_midpoint"] = sorted_df["peaks_midpoint"].round().astype(int)
    sorted_df["start"] = sorted_df["peaks_midpoint"]
    sorted_df["end"] = (sorted_df["peaks_midpoint"] + 1)
    sorted_df = sorted_df.drop_duplicates(["chrom", "start", "end"]) # if any
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


""" For fraction bind with all expressed-and-nonexpressed genes"""

def main():
    """
        Function main to call sub-functions and automate all
    """
    tfs_filelist = glob(tf_files) # file list from TF filepath regex
    tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 3000, 3000)
    final_output_file = "TF_fraction_bind_to_gencodeTSS_3kb_merged.txt"
    
    # Merging of tpm datasets
    tss3kb_pybed = pybedtools.BedTool.from_dataframe(tss_coord_df)
    merged_tss3kb_pybed = tss3kb_pybed.merge(c=[7,7], o=["count_distinct","collapse"])
    merged_tss3kb_df = pd.read_csv(merged_tss3kb_pybed.fn, sep="\t", header=None)
    merged_tss3kb_df.columns = ["chrom", "chrom_Start", "chrom_End", "distinct_count", "gene_symbol"]

    plot_dataset = generate_peaks_singlebinned_perc_bind(merged_tss3kb_df, tfs_filelist)

    return(plot_dataset)


""" For promoter and non-promoter split of peaks """

def main():

    tfs_filelist = glob(tf_files) # file list from TF filepath regex
    tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 3000, 1000)

    # Read promoter (-3kB/+1KB), and filter promoter for non-promoters:
    tss_pybed = pybedtools.BedTool.from_dataframe(tss_coord_df)
    
    if not os.path.exists(join(output_dir, "promoter_based_peaks")):
        os.makedirs(join(output_dir, "promoter_based_peaks"))

    if not os.path.exists(join(output_dir, "nonpromoter_based_peaks")):
        os.makedirs(join(output_dir, "nonpromoter_based_peaks"))

    for each_file in tfs_filelist:
        filename = basename(each_file)
        tf_bedfile = pybedtools.BedTool(each_file)
        print("\nProcessing {}".format(filename))

        # Name output files:
        pybed_outfile = join(join(output_dir, "promoter_based_peaks"), (filename + "_promoter"))
        pybed_outfile_v = join(join(output_dir, "nonpromoter_based_peaks"), (filename + "_nonpromoter"))

        tss_tf_intersect = tf_bedfile.intersect(tss_pybed, u=True, output=pybed_outfile)          
        tss_tf_outersect = tf_bedfile.intersect(tss_pybed, v = True, output=pybed_outfile_v)

if __name__ == '__main__':
    main()


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


def generate_peaks_singlebinned_perc_bind_formerged_tss(tssbins_coord_df, tfs_file_list, **kwargs):

    #tssbins_coord_df = merged_tss3kb_df
    #tfs_file_list = tfs_filelist
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

# For fraction bind with expressed genes:
def main():
    """
        Function main to call sub-functions and automate all
    """
    tfs_filelist = glob(tf_files) # file list from TF filepath regex
    tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 3000, 3000)
    
    # Gene expression data with 1TPM or more:
    gene_exp_df = final_gene_expression_model(join(download_dir,"ENCFF139ZPW.tsv"))
    gene_exp_df = gene_exp_df.loc[:,["gene_id", "TPM", "FPKM"]]
    merged_df = pd.merge(tss_coord_df, gene_exp_df, on="gene_id")
    final_tss_coord_df = merged_df.loc[merged_df["TPM"]>1]
    final_output_file = "TF_fraction_bind_to_gencodeTSS_tpm1_merged.txt"

    # Merging of tpm datasets
    tss3kb_pybed = pybedtools.BedTool.from_dataframe(final_tss_coord_df)
    merged_tss3kb_pybed = tss3kb_pybed.merge(c=[7,7], o=["count_distinct","collapse"])
    merged_tss3kb_df = pd.read_csv(merged_tss3kb_pybed.fn, sep="\t", header=None)
    merged_tss3kb_df.columns = ["chrom", "chrom_Start", "chrom_End", "distinct_count", "gene_symbol"]

    # Fraction bind at TSS:
    plot_dataset = generate_peaks_singlebinned_perc_bind_formerged_tss(merged_tss3kb_df, tfs_filelist)

    return(plot_dataset)

if __name__ == '__main__':
    fraction_bind_data = main()


################################################
# How many TSS falls within 3kb of other TSS: ##
################################################

# For Chris : Merge to find discrete tss promoter:
tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 3000, 3000)
select_cols = ["chrom", "chrom_Start", "chrom_End", "start", "end", "strand", "gene_name"]
tss_coord_df = tss_coord_df.loc[:,select_cols]
tss3kb_pybed = pybedtools.BedTool.from_dataframe(tss_coord_df)
merged_tss3kb_pybed = tss3kb_pybed.merge(c=[7,7], o=["count_distinct","collapse"])
merged_tss3kb_df = pd.read_csv(merged_tss3kb_pybed.fn, sep="\t", header=None)
merged_tss3kb_df.to_csv(join(output_dir, "3kb_TSS_genecode_merged_discrete_promoters.bed"), sep="\t", header=True, index=False)

# Concat TFs file
tfs_filelist = glob(tf_files)
concat_list = []
for each_file in tfs_filelist:
    tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
    print("processing : {}\n".format(tf_name))
    tf_df = pd.read_csv(each_file, sep="\t", header=None)
    tf_df["tf_name"] = tf_name
    concat_list.append(tf_df)

combined_tf_df = pd.concat(concat_list, ignore_index=True)
combined_tf_df.to_csv(join(output_dir,"combined_tf_df.bed"), sep="\t", header=None, index=False)

# download HepG2 FPKM file from encode on fly:
download_dir = output_dir
if not os.path.exists(join(download_dir,"ENCFF139ZPW.tsv")):
    os.system("wget -P {} https://www.encodeproject.org/files/ENCFF139ZPW/@@download/ENCFF139ZPW.tsv".format(download_dir))

def final_gene_expression_model(gene_expression_file):
    expression_df = pd.read_csv(gene_expression_file, sep="\t")
    expression_df["geneID_ensemble"] = expression_df["gene_id"].str.replace(r"\..*", ""); expression_df.head()
    return(expression_df)


def generate_singlebinned_tss_coords_genemodel_with_id(genecode_edited_file, upstream_range, downstream_range):
    """ 
    Function to filter raw genecode model for strand specific coordinate settings 
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
    gencode_pos_df = df_pos.loc[:, select_cols]
    gencode_neg_df = df_neg.loc[:, select_cols]
    gencode_model_df = pd.concat([gencode_pos_df, gencode_neg_df]).sort_values(["chrom", "chrom_Start","chrom_End"])

    # Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    gencode_model_df.loc[gencode_model_df["chrom_Start"] < 0, "chrom_Start"] = 1 # assign start_coords = 1
    gencode_model_df["geneID_gencode"] = gencode_model_df["gene_id"].str.replace(r"\..*", "")
    # tss_coord_df = genecode_model_df.loc[genecode_model_df["chrom_Start"] > 0, :]
    # tss_coord_df.to_csv(join(output_dir, "genecodeV19_singlebinnned_tss_coordinate_info.bed"), sep="\t", index=False, header=False)
    
    return(gencode_model_df)

gene_exp_df = final_gene_expression_model(join(download_dir,"ENCFF139ZPW.tsv")); gene_exp_df.head(700).tail(5)
gencode_tss_coord_df = generate_singlebinned_tss_coords_genemodel_with_id(genecode_file_edited, 3000, 3000)

merged_tss_genexp = pd.merge(gencode_tss_coord_df, gene_exp_df[["geneID_ensemble", "TPM", "FPKM"]], left_on = "geneID_gencode", right_on="geneID_ensemble", how="left")
merged_tss_geneexp_1tpm = merged_tss_genexp.loc[merged_tss_genexp["TPM"] >= 1 ]
select_cols = ["chrom", "chrom_Start", "chrom_End", "start", "end", "strand", "gene_name"]
merged_tss_geneexp_1tpm_df = merged_tss_geneexp_1tpm.loc[:,select_cols]
merged_tss_geneexp_1tpm_df.to_csv(join(output_dir,"final_tss_geneexp_1tpm_df.bed"), sep="\t", header=None, index=False)

merged_tss_geneexp_1tpm_pybed = pybedtools.BedTool.from_dataframe(merged_tss_geneexp_1tpm_df)
final_merged_pybed = merged_tss_geneexp_1tpm_pybed.merge(c=[7,7], o=["count_distinct","collapse"])
final_merged_pybed_df = pd.read_csv(final_merged_pybed.fn, sep="\t", header=None)
final_merged_pybed_df.to_csv(join(output_dir,"final_tss_geneexp_1tpm_df_merged.bed"), sep="\t", header=None, index=False)
