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
final_output_file = "all_tf_tss_intersect.bed"

# if not os.path.exists(output_dir):
#	os.makedirs(output_dir)

# Useful for job submission - if multiple threads running with a race condition to create the dir:
try:
	os.makedirs(suboutput_dir)
except OSError as exc:
	if exc.errno != errno.EEXIST:
		raise
	pass


def final_genecode_model(genecode_edited_file):
	""" 
	Function to filter raw genecode model for strand specific coordinate settings 
	
	Args:
		data_raw(TSV bedfile): Raw genecodeV19 gtf formatted bed file with selective cols
	
	Returns:
		data(pd.DataFame): Processed and cleansed genecode model data
	"""

	# Read Tab-separated value bedfiles:
	df = pd.read_csv(genecode_file_edited, sep="\t")
	regex_pat = re.compile(r";", flags=0) # flags=re.IGNORECASE
	df["gene_name"] = df["gene_name"].str.replace(regex_pat, "")
	df["gene_id"] = df["gene_id"].str.replace(regex_pat, "")
	df_pos = df.loc[df["strand"] == "+"]; df_pos.head()
	df_neg = df.loc[df["strand"] == "-"]; df_neg.head()
	df_neg_new = df_neg.loc[ :,["chrom","end","start", "gene_name", \
							"gene_id", "strand", "feature"]]; df_neg_new.head()
	df_neg_new.columns = df_pos.columns
	df_final = pd.concat([df_pos, df_neg_new]).sort_values(["chrom", "start","end"])
	print("Current dimension of the genecode model:{}".format(df_final.shape))
	final_genecode_model = df_final.drop_duplicates()
	print("Dropping duplicates,if any - current dimension of genecode model:{}\n\n".format(df_final.shape))
	final_genecode_model.to_csv(join(output_dir, "GenecodeV19_gtf_model_final.bed"), sep="\t", header = True, index= False)
	
	return(final_genecode_model)


def generate_tss_binned_coords(genecode_coordinates_info, upstream_range, downstream_range, bin_size):

	""" 
	Function to generate equal bins up-downstream of Tss coordinates from given gene model

	Args:
		genecode_coordinates_info(pd.DataFrame): Pandas Dataframe made for final genomic coords 
		upstream_range(int): Upstream range of Tss for which analysis is requested
		downstream_range(int): Downstream range of Tss for which analysis is requested		
		bin_size(int): Window or bin size requested for analysis

	Returns:
		data(pd.DataFrame): With equally spaced binned coordinates up-downstream of Tsss
	"""

	# Sorting coordinates for consistency required later on: 
	genecode_df =  genecode_coordinates_info.sort_values(["chrom","start","end"])
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  genecode_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	genecode_concat_df = pd.concat([genecode_df]*bin_len, ignore_index="TRUE")
	genecode_sorted_df = genecode_concat_df.sort_values(["chrom","start","end"])

	# Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	# Combine the genecode df and bin df by cbind or column wise:
	temp_genecode_df = pd.concat([genecode_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_genecode_df["tss_midpoint"] = temp_genecode_df["start"]
	select_cols = ["chrom", "tss_midpoint", "bin_start","bin_end", "gene_name","strand", "gene_id", "feature"]
	final_genecode_df = temp_genecode_df.loc[:,select_cols]

	""" For positive strand, i.e if strand1 == "+" 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	genecode_pos_df = final_genecode_df.loc[final_genecode_df["strand"] == "+",]
	genecode_pos_df["chrom_start"] = genecode_pos_df["tss_midpoint"] + genecode_pos_df["bin_start"]
	genecode_pos_df["chrom_end"] = genecode_pos_df["tss_midpoint"] + genecode_pos_df["bin_end"]

	""" For negative strand, i.e if strand1 == "-" 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain higher coords for end sites """

	genecode_neg_df = final_genecode_df.loc[final_genecode_df["strand"] == "-",]
	genecode_neg_df["chrom_start"] = genecode_neg_df["tss_midpoint"] - genecode_neg_df["bin_end"]
	genecode_neg_df["chrom_end"] = genecode_neg_df["tss_midpoint"] - genecode_neg_df["bin_start"]

	# Combine positive and negative stranded genes:
	genecode_model_df = pd.concat([genecode_pos_df, genecode_neg_df])
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'strand', u'gene_name', u'tss_midpoint']
	tss_coord_df = genecode_model_df.loc[:,select_cols]

	# Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
	tss_coord_df = tss_coord_df.loc[tss_coord_df["chrom_start"] > 0, :]
	tss_coord_df.to_csv(join(output_dir, "tss_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(tss_coord_df)


def generate_tss_single_binned_coords(genecode_coordinates_info, upstream_range, downstream_range, bin_size):

	""" 
	Function to generate equal bins up-downstream of Tss coordinates from given gene model

	Args:
		genecode_coordinates_info(pd.DataFrame): Pandas Dataframe made for final genomic coords 
		upstream_range(int): Upstream range of Tss for which analysis is requested
		downstream_range(int): Downstream range of Tss for which analysis is requested		
		bin_size(int): Window or bin size requested for analysis

	Returns:
		data(pd.DataFrame): With equally spaced binned coordinates up-downstream of Tsss
	"""

	# Sorting coordinates for consistency required later on: 
	genecode_df =  genecode_coordinates_info.sort_values(["chrom","start","end"])
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  genecode_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	genecode_concat_df = pd.concat([genecode_df]*bin_len, ignore_index="TRUE")
	genecode_sorted_df = genecode_concat_df.sort_values(["chrom","start","end"])

	# Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	# Combine the genecode df and bin df by cbind or column wise:
	temp_genecode_df = pd.concat([genecode_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_genecode_df["tss_midpoint"] = temp_genecode_df["start"]
	select_cols = ["chrom", "tss_midpoint", "bin_start","bin_end", "gene_name","strand", "gene_id", "feature"]
	final_genecode_df = temp_genecode_df.loc[:,select_cols]

	""" For positive strand, i.e if strand1 == "+" 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	genecode_pos_df = final_genecode_df.loc[final_genecode_df["strand"] == "+",]
	genecode_pos_df["chrom_start"] = genecode_pos_df["tss_midpoint"] + genecode_pos_df["bin_start"]
	genecode_pos_df["chrom_end"] = genecode_pos_df["tss_midpoint"] + genecode_pos_df["bin_end"]

	""" For negative strand, i.e if strand1 == "-" 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain higher coords for end sites """

	genecode_neg_df = final_genecode_df.loc[final_genecode_df["strand"] == "-",]
	genecode_neg_df["chrom_start"] = genecode_neg_df["tss_midpoint"] - genecode_neg_df["bin_end"]
	genecode_neg_df["chrom_end"] = genecode_neg_df["tss_midpoint"] - genecode_neg_df["bin_start"]

	# Combine positive and negative stranded genes:
	genecode_model_df = pd.concat([genecode_pos_df, genecode_neg_df])
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'strand', u'gene_name', u'tss_midpoint']
	tss_coord_df = genecode_model_df.loc[:,select_cols]

	# Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
	tss_coord_df = tss_coord_df.loc[tss_coord_df["chrom_start"] > 0, :]
	tss_coord_df.to_csv(join(output_dir, "tss_coordinate_info.bed"), sep="\t", index=False, header=False)

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


def generate_peaks_binned_perc_bind(tssbins_coord_df, tfs_file_list, **kwargs):
	"""
		Function to generate fraction of peaks binding in each bins up-downstream of Tss
		
		Args:
			tssbins_coord_df(pd.DataFrame): Genecode Tss binned dataframe
			tfs_file_list(pythonlist): List of TF peak files with full path
		
		Returns:
			Fraction of TF binding peaks at each Tss bins requested 
	"""		
	print("kwargs: {}\n".format(kwargs)) # file_name =  kwargs["files_basename"]
	tss_sorted_df = tssbins_coord_df.sort_values(["chrom", "chrom_start", "chrom_end"]).reset_index(drop=True)
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
		# print(tss_tf_intersect.head())  
		# tf_bedfile.intersect(peaks_bedfile, wa = True, wb = True, v = True, output=pybed_outfile_v)

		# Working with the dataframes; reading output file of pybedtool intersect:
		tss_tf_df = pd.read_csv(tss_tf_intersect.fn, sep = "\t", header = None)
		tss_tf_df.rename(columns={3:"bin_start", 4:"bin_end", 8:"peak_chrom", 9:"peak_start", 10:"peak_end"}, inplace=True)
		# tss_tf_df = tss_tf_df.drop_duplicates(["peak_chrom", "peak_start", "end"]).reset_index(drop=True)
		
		# Filter columns based on strings attribute from other numericals:
		select_cols = tss_tf_df.columns[tss_tf_df.columns.map(lambda x: isinstance(x, str)).tolist()]
		final_df = tss_tf_df.loc[:,select_cols]
		final_df["peak_count"] = 1
		df_grouped =  final_df.groupby(["bin_start", "bin_end"]).apply(lambda x : x["peak_count"].sum()/float(tf_bedfile.count()))
		print("Dimension of currently intersected peak file is: {}".format(df_grouped.reset_index().shape))

		if file_basename not in master_dict:
		  master_dict[file_basename] = df_grouped
		print("Intersection of {} completed!!!...\n\n".format(file_basename))

	# Combine all dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	tss_tf_intersect_combined_df = pd.concat(master_dict).reset_index()
	print("Combined dataframe shape : {}".format(tss_tf_intersect_combined_df.shape))
	tss_tf_intersect_combined_df.rename(columns={"level_0":"tf_name", 0:"bind_fraction"}, inplace=True)
	tss_tf_intersect_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)
	
	return(tss_tf_intersect_combined_df)


def z_score(df):
	"""
		Function to compute z-score for value normalization
		
		Args:
			df(pd.DataFrame): Dataframe for which z-score is requested

		Returns:
			pd.Dataframe with z-scores
	"""
	df.columns = [x for x in df.columns.tolist()] # add suffix to x(if needed)
	return ((df - df.mean())/df.std(ddof=0))


def main():
	"""
		Function main to call sub-functions and automate all
		
		Returns:
			Combined heatmap Dataframe with bind fraction for each given TF file, and for each bin  
	"""
	tfs_filelist = glob(tf_files) # file list from TF filepath regex
	genecode_coordinate_df = final_genecode_model(genecode_file_edited)
	tss_coord_df = generate_tss_binned_coords(genecode_coordinate_df, 5000, 5000, 500)
	plot_dataset = generate_peaks_binned_perc_bind(tss_coord_df, tfs_filelist )
	tf_tss_heatplot_df = plot_dataset.pivot(index="tf_name", columns="bin_start", values="bind_fraction")
	tf_tss_heatplot_df.fillna(0, inplace=True) # points to no-binding event
	tf_tss_heatplot_df.to_csv(join(output_dir, "all_tf_tss_intersect_for_heatmap.txt"), sep ="\t", header = True, index = True)
	tf_tss_heatplot_df.to_pickle(join(output_dir, "all_tf_tss_intersect_for_heatmap.pkl")) # to read: pd.read_pickel(filename.pkl)
	
	# Z-score normalized binding fraction heatmap dataframe
	tf_tss_zscore_heatplot_df = z_score(tf_tss_heatplot_df)
	tf_tss_zscore_heatplot_df.to_csv(join(output_dir, "all_tf_tss_intersect_for_heatmap_zscore.txt"), sep ="\t", header = True, index = True)
	print("Time taken to complete analyis: {}".format(time.time()-start_time))

	return(tf_tss_heatplot_df, tf_tss_zscore_heatplot_df)

if __name__ == '__main__':
	final_data, final_data_zscore_norm = main()



######################################################################
# Simple genecode analysis expanding 3kb up-downstream with no bins
######################################################################


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
#	os.makedirs(output_dir)

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

# Top 5 genes with highest TPM vals:
merged_tss_genexp.nlargest(5, ["TPM"])

""" Local Machine computation """
input_dir = '/Users/suryachhetri/Dropbox/for_genemodels/number_for_chris'

#Case 0 
#####################
combined_tf_df_sorted = pd.read_csv(join(input_dir, "combined_tf_df_sorted.bed"), sep="\t", header=None)
combined_tf_df_sorted_pybed = pybedtools.BedTool.from_dataframe(combined_tf_df_sorted)
####################

# Case 1
merged_discrete_tss_df = pd.read_csv(join(input_dir, "3kb_TSS_genecode_merged_discrete_promoters.bed"), sep="\t")
merged_discrete_tss_pybed = pybedtools.BedTool.from_dataframe(merged_discrete_tss_df)

# Case 0/1 intersect
tss_with_tfpeaks = merged_discrete_tss_pybed.intersect(combined_tf_df_sorted_pybed, u=True)
tss_with_0_tfpeaks = merged_discrete_tss_pybed.intersect(combined_tf_df_sorted_pybed, v=True)
(tss_with_tfpeaks.count(), tss_with_0_tfpeaks.count())
# tss_with_0_tfpeaks = pd.read_csv(tss_with_0_tfpeaks.fn, sep="\t", header=None)

# Case 2
final_tss_geneexp_1tpm_df = pd.read_csv(join(input_dir, "final_tss_geneexp_1tpm_df.bed"), sep="\t", header=None)
final_tss_geneexp_1tpm_pybed = pybedtools.BedTool.from_dataframe(final_tss_geneexp_1tpm_df)

# Case 0/2 intersect
tss_with_tfpeaks = final_tss_geneexp_1tpm_pybed.intersect(combined_tf_df_sorted_pybed, u=True)
tss_with_0_tfpeaks = final_tss_geneexp_1tpm_pybed.intersect(combined_tf_df_sorted_pybed, v=True)
(tss_with_tfpeaks.count(), tss_with_0_tfpeaks.count())
# tss_with_0_tfpeaks = pd.read_csv(tss_with_0_tfpeaks.fn, sep="\t", header=None)

# Case 3
final_tss_geneexp_1tpm_df_merged = pd.read_csv(join(input_dir, "final_tss_geneexp_1tpm_df_merged.bed"), sep="\t", header=None)
final_tss_geneexp_1tpm_merged_pybed = pybedtools.BedTool.from_dataframe(final_tss_geneexp_1tpm_df_merged)

# Case 0/3 intersect
tss_with_tfpeaks = final_tss_geneexp_1tpm_merged_pybed.intersect(combined_tf_df_sorted_pybed, u=True)
tss_with_0_tfpeaks = final_tss_geneexp_1tpm_merged_pybed.intersect(combined_tf_df_sorted_pybed, v=True)
(tss_with_tfpeaks.count(), tss_with_0_tfpeaks.count())
# tss_with_0_tfpeaks = pd.read_csv(tss_with_0_tfpeaks.fn, sep="\t", header=None)

################################################################
################################################################

select_cols = ["chrom", "start", "end", "strand", "gene_name"]
tss_df = tss_coord_df.loc[:,select_cols]
df_pos = tss_df.loc[tss_df["strand"] == "+"]; df_pos.head()
df_neg = tss_df.loc[tss_df["strand"] == "-"]; df_neg.head()
# Maintain chrom_start coord > chrom_neg to avoid bed-intersection upset 
df_neg.rename(columns={"start" : "end", "end" : "start"}, inplace=True)

# Combine positive and negative stranded genes:
select_cols = ["chrom", "start", "strand", "gene_name"]
genecode_pos_df = df_pos.loc[:, select_cols]
genecode_neg_df = df_neg.loc[:, select_cols]
genecode_model_df = pd.concat([genecode_pos_df, genecode_neg_df]).sort_values(["chrom", "start"])
genecode_model_df["end"] = genecode_model_df["start"] + 1

# Rearrange columns:
select_cols = ["chrom", "start", "end", "strand", "gene_name"]
tssfinal_df = genecode_model_df.loc[:,select_cols]
tss_pybed = pybedtools.BedTool.from_dataframe(tssfinal_df)

tss_intersect = tss_pybed.intersect(tss3kb_pybed, u=True)
tss_intersect_df = pd.read_csv(tss_intersect.fn, header=None, sep="\t")
tss_intersect_df.head()
tss_intersect_df.drop_duplicates([0,1,2])

# map tss to bed region and count with col 5:
tss_coord_df = generate_singlebinned_tss_coords_genemodel(genecode_file_edited, 3000, 3000)
tss_coord_df["count"] = 1
tss_coord_df["count"] = tss_coord_df["count"].astype(int)
select_cols = ["chrom", "chrom_Start", "chrom_End", "gene_name", "count", "strand"]
tss_coord_df = tss_coord_df.loc[:,select_cols]
tss3kb_pybed = pybedtools.BedTool.from_dataframe(tss_coord_df)

select_cols = ["chrom", "start", "end", "strand", "gene_name"]
tss_df = tss_coord_df.loc[:,select_cols]
df_pos = tss_df.loc[tss_df["strand"] == "+"]; df_pos.head()
df_neg = tss_df.loc[tss_df["strand"] == "-"]; df_neg.head()
# Maintain chrom_start coord > chrom_neg to avoid bed-intersection upset 
df_neg.rename(columns={"start" : "end", "end" : "start"}, inplace=True)

# Combine positive and negative stranded genes:
select_cols = ["chrom", "start", "strand", "gene_name"]
genecode_pos_df = df_pos.loc[:, select_cols]
genecode_neg_df = df_neg.loc[:, select_cols]
genecode_model_df = pd.concat([genecode_pos_df, genecode_neg_df]).sort_values(["chrom", "start"])
genecode_model_df["end"] = genecode_model_df["start"] + 1

genecode_model_df["count"] = 1
genecode_model_df["count"] = genecode_model_df["count"].astype(int)
select_cols = ["chrom", "start", "end", "gene_name", "count", "strand"]
tssfinal_df = genecode_model_df.loc[:,select_cols]
tss_pybed = pybedtools.BedTool.from_dataframe(tssfinal_df)

tss_map = tss3kb_pybed.map(tss_pybed, wa=True, wb=True, c=[5,4], o=["sum","collapse"], null="null")
tss_map_df = pd.read_csv(tss_map.fn, header=None, sep="\t")
tss_map_df.drop_duplicates([0,1,2])



##################################################
# python ggploting locally : (Gencode Tss barplot)
##################################################


# Run @ local computer:
output_dir = ("~/for_chris/batch_I/for_genecodeTssProfile/for_genemodels")
final_output_file = 'all_tf_fraction_bind_to_gencodeTSS.txt'
copy_file = join(output_dir, final_output_file)	
!scp newcluster:$copy_file .

import matplotlib; import matplotlib.pyplot as plt
from ggplot import *
from plotnine import *
import pandas as pd, numpy as np

# matplotlib.get_backend()
# matplotlib.use("Qt4Agg") # remote

df_tss = pd.read_csv('all_tf_fraction_bind_to_gencodeTSS.txt', sep="\t")
df_tss["log2_peak_count"] = np.log2(df_tss["peak_count"])
df_tss.to_csv("all_tf_fraction_binding_gencodeTSS_final.txt", sep ="\t", header = True, index = True)

# Ordering of factor levels:
sorted_df_tss = df_tss.sort_values(["bind_fraction"]).reset_index(drop=True)
df_tss['tf_name'] = pd.Categorical(df_tss['tf_name'], categories=sorted_df_tss["tf_name"], ordered=True)

# fifty_perc = sorted_df_tss.loc[sorted_df_tss["bind_fraction"]>0.5]
# fifty_perc.to_csv("all_tf_fraction_binding_gencodeTSS_50percent_final.txt", sep ="\t", header = True, index = True)

plot = (ggplot(df_tss) + 
		aes("tf_name", "bind_fraction", color="log2_peak_count" ) + 
		geom_bar(stat="identity") +
		coord_flip() +
		ggtitle('Binding fraction of 208 factors at Gencode TSS (n=57820)') +
		scale_x_discrete(name="Transcription Factors") +  scale_y_continuous(name="Binding Fraction") +
		# xlab("Transcription Factors") +  ylab("Binding Fraction") +
		guides(fill=guide_legend(title="Log2(Peak Count)")) +
		geom_hline(aes(yintercept=0.5), size=0.55, color="black", linetype="dashed") + 
		theme_minimal() + 
		theme(
		axis_text_y = element_text(size=1.3), #theme(axis_text_x=element_text(angle=45))
		plot_title = element_text(size=9, face="bold", hjust = 0.6),
		legend_title = element_text(size=8, face="bold") 
		) 
	)

ggsave(plot, "TFbind_tss_3kb_genecode_barplot_fig.pdf", width = 6, height = 5.5) 



#########################################################################
# IDEAS piechart analysis-discrepancy impact-ratio and heatmap clustering
#########################################################################


import pandas as pd, numpy as np
from glob import glob
import re
from os.path import join, basename, dirname

file_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/piecharts_ideas_unique_newHepG2/*intersectbed_counts_with_ideas.txt")
df_list = []
tf_list = []

regex_pat = re.compile(r"_intersectbed_counts_with_ideas.txt")
for each_file in file_list:
	tf_name = regex_pat.split(basename(each_file))[0]
	df = pd.read_csv(each_file, sep="\t", header=None)
	df.columns = ["state", "peak_count"]
	df["tf_name"] = tf_name
	df = df.loc[:,["tf_name", "state", "peak_count"]]
	df_list.append(df)
	tf_list.append(tf_name)

combined_df = pd.concat(df_list, ignore_index=True)
final_df = combined_df.sort_values(["tf_name"]).reset_index(drop=True)
sum_df = final_df.groupby(["tf_name"])["peak_count"].sum().reset_index(drop=False)
final_df["total_count"] = final_df.groupby(["tf_name"])["peak_count"].transform(sum)
final_df["fraction_bind"] = final_df["peak_count"].astype(float)/final_df["total_count"]
heatmap_df = final_df.pivot(index="tf_name", columns="state", values="fraction_bind")
heatmap_df.to_csv(join(output_dir, "ideas_piestate_dist_for_heatmap.txt"), sep ="\t", header = True, index = True)
# count_df = final_df.groupby(["tf_name"])["total_count"].value_counts()

peak_tf_list = []
peak_count_list = []
peak_file_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*narrowPeak*")
regex_pat = re.compile(r"narrowPeak_")
for each_file in peak_file_list:
	tf_name = regex_pat.split(basename(each_file))[1]
	df = pd.read_csv(each_file, sep="\t", header=None)
	peak_count_list.append(df.shape[0])
	peak_tf_list.append(tf_name)

peak_df = pd.DataFrame({"tf_name": peak_tf_list, 
						"peak_counts" : peak_count_list })
peak_df = peak_df.loc[:,["tf_name","peak_counts"]]
peak_df = peak_df.sort_values(["tf_name"]).reset_index(drop=True)
peak_df.to_csv(join(output_dir, "all_tfs_peakcount.txt"), sep ="\t", header = True, index = False)

df = pd.concat([peak_df, sum_df], axis=1)
df["diff"] = df["peak_counts"] - df["peak_count"]
df["impact_ratio"] = df["diff"].astype(float)/df["peak_counts"].astype(float)

# Explaination note: All the peaks with less than 50% overlap is disqualified for peak distribution
diff_df = df.iloc[df["diff"].nlargest(10).index]

# The impact ratio of disqualified peaks shows that max impact would be of CBX1, with shift of just 0.4%(0.004831*100) 
# at max, that is percentage shift it might have caused - had those disqualified peaks been included in the piechart. 
impact_df = df.iloc[df["impact_ratio"].nlargest(10).index]


######################################################################################################
# IDEAS piechart data (Data processing for stacked barplot that only has promoter-enhancer assoc dist:
######################################################################################################


import pandas as pd, numpy as np
from glob import glob
import re, os
from os.path import join, basename, dirname

output_dir = os.path.expanduser("~/for_chris/batch_I/for_genecodeTssProfile/for_genemodels")
file_list = glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/piecharts_ideas_unique_newHepG2/*intersectbed_counts_with_ideas.txt")

df_list = []
tf_list = []
regex_pat = re.compile(r"_intersectbed_counts_with_ideas.txt")
for each_file in file_list:
	tf_name = regex_pat.split(basename(each_file))[0]
	df = pd.read_csv(each_file, sep="\t", header=None)
	df.columns = ["state", "peak_count"]
	df["tf_name"] = tf_name
	df = df.loc[:,["tf_name", "state", "peak_count"]]
	df_list.append(df)
	tf_list.append(tf_name)

combined_df = pd.concat(df_list, ignore_index=True)
final_df = combined_df.sort_values(["tf_name"]).reset_index(drop=True)

### For whole genome annotation:
map_ideas = {"Tss" : "Prom assoc", "TssF" :  "Prom assoc", 'TssCtcf': "Prom assoc", 'PromCtcf': "Prom assoc", 
"TssW" : "Prom assoc", "PromP": "Prom assoc", 'PromF1': "Prom assoc", 'PromF2': "Prom assoc",
"Pol2" : "Gene body", "Gen5" : "Gene body", "Gen3": "Gene body", "Gen3Ctcf" : "Gene body", "Elon" : "Gene body", "ElonW" : "Gene body",
"CtcfO" : "Ctcf assoc", "Ctcf" : "Ctcf assoc", 
"Enh" : "Strong Enh", "EnhF" : "Strong Enh", "EnhWF1" : "Weak Enh", "EnhWF2" : "Weak Enh", "EnhWF3" : "Weak Enh", "EnhW" : "Weak Enh",
"Repr1" : "Heterochrom_repressed", "Repr2" : "Heterochrom_repressed", "ReprD":  "Heterochrom_repressed", "LowReprW":  "Heterochrom_repressed", 
"Low1":  "Heterochrom_repressed", "Low2":  "Heterochrom_repressed", "Quies" :  "Heterochrom_repressed", "Art1":  "Heterochrom_repressed",
"Art2":  "Heterochrom_repressed", "Zero" : "Heterochrom_repressed", 
"DnaseD1": "Euchromatin assoc", "DnaseD2": "Euchromatin assoc", "FaireW1": "Euchromatin assoc", "FaireW2": "Euchromatin assoc"}

final_df["ideas_state"] =  final_df["state"].map(map_ideas)
final_df["total_count"] = final_df.groupby(["tf_name"])["peak_count"].transform(sum)
prom_enh_list = ["Prom assoc", "Strong Enh", "Weak Enh"]
prom_enh_df = final_df.loc[final_df["ideas_state"].isin(prom_enh_list)]

remap_ideas = {"Prom assoc" : "Promoter associated", "Weak Enh" : "Enhancer associated", "Strong Enh" : "Enhancer associated"}
prom_enh_df["ideas_state_new"] = prom_enh_df["ideas_state"].map(remap_ideas)
prom_enh_grouped = prom_enh_df.groupby(["tf_name", "ideas_state_new","total_count"])["peak_count"].apply(lambda x: x.sum())
prom_enh_grouped = prom_enh_grouped.reset_index()
prom_enh_grouped["bind_fraction"] = prom_enh_grouped["peak_count"]/prom_enh_grouped["total_count"]
# fraction_sum_test = prom_enh_grouped.groupby(["tf_name"])["bind_fraction"].apply(lambda x: x.sum()).reset_index()
prom_enh_grouped.to_csv(join(output_dir, "ideas_stacked_prom_enh_piestate_dist_barplot_data.txt"), sep ="\t", header = True, index = False)


######################################################################################
# Python ggploting locally : "ideas_stacked_prom_enh_piestate_dist_barplot_data.txt"
######################################################################################

import matplotlib; import matplotlib.pyplot as plt
from ggplot import *
from plotnine import *
import pandas as pd, numpy as np

# matplotlib.get_backend()
# matplotlib.use("Qt4Agg") # remote

df_cis = pd.read_csv('ideas_stacked_prom_enh_piestate_dist_barplot_data.txt', sep="\t")
df_cis = df_cis.loc[:,["tf_name", "ideas_state_new", "peak_count", "total_count", "bind_fraction"]]
df_cis["log2_peak_count"] = np.log2(df_cis["peak_count"])

# Align TFs to previous TF list(made for TSS barplot):
df_tss = pd.read_csv('all_tf_fraction_bind_to_gencodeTSS.txt', sep="\t")
df_tss["log2_peak_count"] = np.log2(df_tss["peak_count"])

# Ordering of factor levels wrt to previous tss based TF list:
sorted_df_tss = df_tss.sort_values(["bind_fraction"]).reset_index(drop=True)
df_cis['tf_name'] = pd.Categorical(df_cis['tf_name'], categories=sorted_df_tss["tf_name"], ordered=True)

# Order the factor/categorical variable to color legend accordingly:
df_cis["ideas_state_new"] = pd.Categorical(df_cis["ideas_state_new"], categories=["Promoter associated", "Enhancer associated"], ordered=True)

# Reorder the factor/categorical variable to color legend accordingly:
df_cis_new = df_cis[df_cis["ideas_state_new"] == "Promoter associated"]
sorted_df_cis = df_cis_new.sort_values(["bind_fraction"]).reset_index(drop=True)
df_cis['tf_name'] = pd.Categorical(df_cis['tf_name'], categories=sorted_df_cis["tf_name"], ordered=True)
df_cis["ideas_state_new"] = pd.Categorical(df_cis["ideas_state_new"], categories=["Promoter associated", "Enhancer associated"], ordered=True)


# prom_enhancer stacked/facet barplot with custom color:
plot = (ggplot(df_cis) + 
		aes("tf_name", "bind_fraction", fill="ideas_state_new") + 
		geom_bar(stat="identity") +
		coord_flip() +
		ggtitle('Binding fraction of 208 factors across Cis-regulatory regions') +
		scale_x_discrete(name="Transcription Factors") +  scale_y_continuous(name="Binding Fraction") +
		scale_fill_manual(name="IDEAS Anno",values=["red","orange"]) +
		guides(fill=guide_legend(title="Annotation")) +
		geom_hline(aes(yintercept=0.5), size=0.55, color="black", linetype="dashed") + 
		theme_minimal() + 
		theme(
		axis_text_y = element_text(size=1.3), #theme(axis_text_x=element_text(angle=45))
		plot_title = element_text(size=9, face="bold", hjust = 0.6),
		legend_title = element_text(size=8, face="bold") 
		) +
		facet_wrap('~ideas_state_new')
	)

ggsave(plot, "TFbind_prom_enh_facetstack_barplot_fig_redo.pdf", width = 6, height = 5.5) 

# prom_enhancer facet_wrap barplot:
plot = (ggplot(df_cis) + 
		aes("tf_name", "bind_fraction", fill="log2_peak_count") + 
		geom_bar(stat="identity") +
		coord_flip() +
		ggtitle('Binding fraction of 208 factors across Cis-regulatory regions') +
		scale_x_discrete(name="Transcription Factors") +  scale_y_continuous(name="Binding Fraction") +
		# scale_fill_manual(name="IDEAS Anno",values=["red","orange"]) +
		# guides(fill=guide_legend(title="Annotation")) +
		geom_hline(aes(yintercept=0.5), size=0.55, color="black", linetype="dashed") + 
		theme_minimal() + 
		theme(
		axis_text_y = element_text(size=1.3), #theme(axis_text_x=element_text(angle=45))
		plot_title = element_text(size=9, face="bold", hjust = 0.6),
		legend_title = element_text(size=8, face="bold") 
		) +
		facet_wrap('~ideas_state_new')

	)

ggsave(plot, "TFbind_prom_enh_facet_barplot_fig_redo.pdf", width = 6, height = 5.5) 



