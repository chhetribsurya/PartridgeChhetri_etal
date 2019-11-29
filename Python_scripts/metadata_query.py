#!/usr/bin/env python

import pandas as pd
from collections import Counter
from os import makedirs, rmdir, remove
from os.path import join, splitext, basename
from os.path import expanduser, exists
from glob import glob
import requests, json
HEADERS = {'accept': 'application/json'}


input_dir = expanduser("~/Dropbox/for_encode/encode3_paper")
output_dir = expanduser("~/Dropbox/for_encode/encode3_paper/chipseq_metadata_info")
se_liblist_outdir = join(output_dir, "SE_libraries_info")
pe_liblist_outdir = join(output_dir, "PE_libraries_info")
#input_dir = expanduser("~/k562_chip")
#output_dir = expanduser("~/k562_chip")
#se_liblist_outdir = join(output_dir, "SE_libraries_info")
#pe_liblist_outdir = join(output_dir, "PE_libraries_info")

if not exists(output_dir):
	makedirs(output_dir)

if not exists(se_liblist_outdir):
	makedirs(se_liblist_outdir)

if not exists(pe_liblist_outdir):
	makedirs(pe_liblist_outdir)

def metadata_idr_query(metadata_tsv_file, cell_type):
	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	#print "\nHere's the list of features: %s \n\n" %(df.columns)
	final_df = df[(df.loc[:,"Assembly"] == "hg19") & (df.loc[:, "Output type"] == "optimal idr thresholded peaks")]
 	#final_df.loc[:, select_cols].head()
 	select_cols = ["File accession", "File format", "Output type", "Experiment accession", 
 					"Experiment target", "Biosample term name", "Run type", "md5sum", "File download URL"]
	#print "However, only picking up these features:%s\n\n" %(select_cols)
	selected_final_df = final_df.loc[:, select_cols]
 	#selected_final_df.loc[:, select_cols].head()
 	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
 	celltype_df_eGFP = celltype_df[celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
 	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
	eGFP_count = celltype_df_eGFP.loc[:,"Experiment target"].value_counts().shape
 	wo_eGFP_count = celltype_df_wo_eGFP.loc[:,"Experiment target"].value_counts().shape
	total_count = celltype_df.loc[:,"Experiment target"].value_counts().shape
	celltype_df.to_csv(join(output_dir, cell_type + "_final_IDR_metadata.txt"), sep="\t")
	print_data = "\n\n%s final data -\nUnique_factor_count : %s\nWithout_eGFP_count : %s\neGFP_count : %s \n\n" \
				%(cell_type, total_count[0], wo_eGFP_count[0], eGFP_count[0])
	celltype_df["File download URL"].to_csv(join(output_dir, cell_type + "_IDR_download_urls.txt"), index=False)
	print(print_data)

	### Create dictionary - to extract lab name using "Experiment accesion":
	fastq_df = df[df.loc[ :, "File format"] == "fastq"][["Experiment accession", "Lab", "File accession", "File format", "Run type"]]
	fastq_lab_dict = fastq_df.set_index("Experiment accession").to_dict(orient="index")
	unique_non_eGFP_df = celltype_df_wo_eGFP.drop_duplicates(["Experiment target"])
	
	Lab_list = []
	for each_acc in unique_non_eGFP_df["Experiment accession"].tolist():
		lab_id = fastq_lab_dict.get(each_acc).get("Lab")
		Lab_list.append(lab_id)
	non_eGFP_Lab_counts = pd.Series(Lab_list).value_counts()

	### Write comments in CSV file with pandas - open in "a" append mode:
	if exists(join(output_dir, cell_type + "_uniqueTF_Count_withIDR_and_Labinfo.txt")):
		remove(join(output_dir, cell_type + "_uniqueTF_Count_withIDR_and_Labinfo.txt")) # avoids append mode appending each time
	with open(join(output_dir, cell_type + "_uniqueTF_Count_withIDR_and_Labinfo.txt"), "a") as f:
		f.write(print_data)
		#f.write("# Unique TFCount distribution with IDR peaks\n")
		f.write("# Lab distribution\n")
		non_eGFP_Lab_counts.to_csv(f, sep="\t")
	return(non_eGFP_Lab_counts)

idr_metadata_k562 =  metadata_idr_query("chipseq_final_metadata.tsv", "K562")
idr_metadata_GM12878 = metadata_idr_query("chipseq_final_metadata.tsv", "GM12878")
# idr_metadata_HepG2 = metadata_idr_query("chipseq_final_metadata.tsv", "HepG2")
# idr_metadata_HEK293 = metadata_idr_query("chipseq_final_metadata.tsv", "HEK293")
# idr_metadata_A549 = metadata_idr_query("chipseq_final_metadata.tsv", "A549")
# idr_metadata_MCF7 = metadata_idr_query("chipseq_final_metadata.tsv", "MCF-7")
# idr_metadata_HeLaS3 = metadata_idr_query("chipseq_final_metadata.tsv", "HeLa-S3")
# idr_metadata_SKNSH = metadata_idr_query("chipseq_final_metadata.tsv", "SK-N-SH")



def metadata_fastq_query(metadata_tsv_file, cell_type):
	#metadata_tsv_file = "Blood_chipseq_metadata_final.tsv"
	#cell_type = "K562"
	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	#print "\nHere's the list of features: %s \n\n" %(df.columns)
	final_df = df[(df.loc[:, "File format"] == "fastq") & (df.loc[:, "File Status"] != "archived")]
	select_cols = ["File accession", "File format", "Output type", "Experiment accession",
						"Experiment target", "Biosample term name", "Run type", "Lab","md5sum", "File download URL"]
	#print "However, only picking up these features:%s\n\n" %(select_cols)
	selected_final_df = final_df.loc[:, select_cols]
	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
	celltype_df_eGFP = celltype_df[celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
	celltype_df_SE = celltype_df[celltype_df.loc[:,"Run type"].str.contains("single-ended")]
	celltype_df_PE = celltype_df[celltype_df.loc[:,"Run type"].str.contains("paired-ended")]

	## Total experiment or expt. accession count data
	total_count = celltype_df.loc[:,"Experiment accession"].value_counts().shape
	eGFP_count = celltype_df_eGFP.loc[:,"Experiment accession"].value_counts().shape
	wo_eGFP_count = celltype_df_wo_eGFP.loc[:,"Experiment accession"].value_counts().shape
	total_df = celltype_df.drop_duplicates(["Experiment accession"])
	lab_count = total_df.loc[:,"Lab"].value_counts()
	SE_count = celltype_df_SE.loc[:,"Experiment accession"].value_counts().shape
	PE_count = celltype_df_PE.loc[:,"Experiment accession"].value_counts().shape

	#celltype_df.to_csv(join(output_dir, cell_type + "_final_fastq_allTF_metadata.txt"), sep="\t")
	print_acc_data = "\n\n%s final data -\nTotal_factor-experiment_count : %s\nWithout_eGFP_count : %s\neGFP_count : %s\nSE_count : %s\nPE_count : %s\n" \
				%(cell_type, total_count[0], wo_eGFP_count[0], eGFP_count[0], SE_count[0], PE_count[0])
	#celltype_df["File download URL"].to_csv(join(output_dir, cell_type + "_allTF_download_urls.txt"), index=False)
	print(print_acc_data)
	print(lab_count)

	## unique TF count data
	total_count = celltype_df.loc[:,"Experiment target"].drop_duplicates().shape
	eGFP_count = celltype_df_eGFP.loc[:,"Experiment target"].drop_duplicates().shape
	wo_eGFP_count = celltype_df_wo_eGFP.loc[:,"Experiment target"].drop_duplicates().shape
	uniq_df = celltype_df.drop_duplicates(["Experiment target"]).reset_index(drop=True)
	SE_count = uniq_df.loc[uniq_df.loc[:, "Run type"]=="single-ended"].shape
	PE_count = uniq_df.loc[uniq_df.loc[:, "Run type"]=="paired-ended"].shape
	#celltype_df.to_csv(join(output_dir, cell_type + "_final_fastq_uniqueTF_metadata.txt"), sep="\t")
	print_tf_data = "\n\nUnique_factor_count : %s\nWithout_eGFP_count : %s\neGFP_count : %s\nSE_count : %s\nPE_count : %s\n" \
				%(total_count[0], wo_eGFP_count[0], eGFP_count[0], SE_count[0], PE_count[0])
	#celltype_df["File download URL"].to_csv(join(output_dir, cell_type + "_uniqueTF_download_urls.txt"), index=False)
	print(print_tf_data)

	### Write comments in CSV file with pandas - open in "a" append mode:
	if exists(join(output_dir, cell_type + "_final_fastq_allTF_metadata.txt")):
		remove(join(output_dir, cell_type + "_final_fastq_allTF_metadata.txt")) # avoids append mode appending each time
	with open(join(output_dir, cell_type + "_final_fastq_allTF_metadata.txt"), "a") as f:
		f.write(print_acc_data)
		f.write("\n# Lab distribution\n")
		lab_count.to_csv(f, sep="\t")
		f.write(print_tf_data)	
	return()


fastq_metadata_k562 = metadata_fastq_query("Blood_chipseq_metadata_final.tsv", "K562")
fastq_metadata_GM12878 = metadata_fastq_query("Blood_chipseq_metadata_final.tsv", "GM12878")



def create_fastq_urlList_to_download(metadata_tsv_file, cell_type):
	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	#print "\nHere's the list of features: %s \n\n" %(df.columns)
	final_df = df[(df.loc[:, "File format"] == "fastq") & (df.loc[:, "File Status"] != "archived")]
	select_cols = ["File accession", 'Controlled by', "Experiment accession",
						"Experiment target", "Biosample term name"]
	selected_final_df = final_df.loc[:, select_cols]
	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]

	reps_fastq_list = celltype_df_wo_eGFP['File accession'].unique().tolist()
	control_fastq_series = celltype_df_wo_eGFP['Controlled by'].str.replace("/files/", "").str.replace("/", "").apply(lambda x: x.split(",")[0])
	control_fastq_list = control_fastq_series.unique().tolist()
	s1 = pd.Series(reps_fastq_list)
	s2 = pd.Series(control_fastq_list)
	final_list = pd.concat([s1,s2], ignore_index=True).tolist()

	#download_list = [ each for each in final_list if each not in downloaded_file_list ]
	full_url_list = [ ] 
	for each in final_list:
		full_url_list.append(join("https://www.encodeproject.org/files", each, "@@download", each+".fastq.gz"))

	"""
	Transfer this list to morgan for downloading of rest of fastqs
	"""
	download_url_df = pd.Series(full_url_list)
	download_url_df.to_csv(join(output_dir, cell_type+"_all_url_list_for_fastq_download.txt"), sep="\t", header=None, index=False)
	
	return(download_url_df)

fastq_list_k562 = create_fastq_urlList_to_download("Blood_chipseq_metadata_final.tsv", "K562")
fastq_list_GM12878 = create_fastq_urlList_to_download("Blood_chipseq_metadata_final.tsv", "GM12878")



def filter_fastq_from_preexisting_urlList_to_download(metadata_tsv_file, pre_fastq_download_file, cell_type):
	#pre_fastq_download_file="k562_encode_fastq_earlier_downloads_list.txt"
	with open(join(input_dir, pre_fastq_download_file), "r") as f:
		file_list = [basename(url).split(".")[0] for url in f ]

	downloaded_file_list = file_list

	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	final_df = df[(df.loc[:, "File format"] == "fastq") & (df.loc[:, "File Status"] != "archived")]
	select_cols = ["File accession", 'Controlled by', "Experiment accession",
						"Experiment target", "Biosample term name"]
	selected_final_df = final_df.loc[:, select_cols]
	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]

	reps_fastq_list = celltype_df_wo_eGFP['File accession'].unique().tolist()
	control_fastq_series = celltype_df_wo_eGFP['Controlled by'].str.replace("/files/", "").str.replace("/", "").apply(lambda x: x.split(",")[0])
	control_fastq_list = control_fastq_series.unique().tolist()
	s1 = pd.Series(reps_fastq_list)
	s2 = pd.Series(control_fastq_list)
	final_list = pd.concat([s1,s2], ignore_index=True).tolist()

	download_list = [ each for each in final_list if each not in downloaded_file_list ]
	full_url_list = [ ] 
	for each in download_list:
		full_url_list.append(join("https://www.encodeproject.org/files", each, "@@download", each+".fastq.gz"))

	"""
	Transfer this list to morgan for downloading of rest of fastqs
	"""
	download_url_df = pd.Series(full_url_list)
	download_url_df.to_csv(join(output_dir, cell_type+"_filtered_url_from_preDownloaded_list.txt"), sep="\t", header=None, index=False)
	return(download_url_df)

filtered_fastq_list_k562 = filter_fastq_from_preexisting_urlList_to_download("Blood_chipseq_metadata_final.tsv", "k562_encode_fastq_earlier_downloads_list.txt", "K562")
filtered_fastq_list_GM12878 = filter_fastq_from_preexisting_urlList_to_download("Blood_chipseq_metadata_final.tsv", "k562_encode_fastq_earlier_downloads_list.txt", "GM12878")


#cell_type = "GM12878"
#control_metadata_tsv_file = "control_Blood_chipseq_metadata_final.tsv"
cell_type = "K562"
metadata_tsv_file = "Blood_chipseq_metadata_final.tsv"
fastq_source_lib_dir = join(input_dir, "fastqs")
fastq_dest_lib_dir = join(input_dir, "k562_SE_fastq_libraries")

def structure_SE_fastq_libraries_for_peakcall(metadata_tsv_file, cell_type, fastq_source_lib_dir, fastq_dest_lib_dir):
	if not exists(fastq_dest_lib_dir):
		makedirs(fastq_dest_lib_dir)

	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	""" Filter Ended fastq files """
	#df_pe = df[df.loc[:, "Run type"] == "paired-ended"]
	df_se = df[df.loc[:, "Run type"] == "single-ended"]
	final_df = df_se[(df_se.loc[:, "File format"] == "fastq") & (df_se.loc[:, "File Status"] != "archived")]
	select_cols = ["File accession", 'Controlled by', "Experiment accession", "Biological replicate(s)", "Lab",
							"Experiment target", "Biosample term name"]
	selected_final_df = final_df.loc[:, select_cols]
	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
	select_cols = ['File accession', "Experiment target", 'Controlled by', 'Paired end','Experiment accession', 'Biological replicate(s)']
	sorted_dframe = celltype_df_wo_eGFP.sort_values(["Experiment accession", "Biological replicate(s)"])
	sorted_dframe = sorted_dframe.reset_index(drop=True)

	### Select the samples with only controls available:
	sorted_dframe = sorted_dframe[~sorted_dframe['Controlled by'].isnull()]
	sorted_dframe["Controlled by"] = sorted_dframe['Controlled by'].replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True) # capture regex

	### Filter out sample with 1 reps or selecting samples with 2 or more reps:
	rep_idx = []
	expt_acc_list = sorted_dframe.loc[:,"Experiment accession"].tolist()
	for i,item in enumerate(expt_acc_list):
		if (expt_acc_list.count(item) >= 2):
			#print(item), i 
			rep_idx.append(i)

	final_rep_df = sorted_dframe.iloc[rep_idx]
	final_rep_df = final_rep_df.reset_index(drop=True)

	### Select top(head-2) replicates from the list (if more than 2 reps for any experiment):
	acc_grouped_df = final_rep_df.groupby(["Experiment accession"]).head(2).reset_index(drop=True)
	acc_grouped_df["Experiment_acc_repId"] = final_rep_df[["Experiment accession", "Biological replicate(s)"]].astype(str).apply(lambda x: '_'.join(x), axis=1)
	acc_grouped_df = acc_grouped_df.sort_values(["Experiment_acc_repId"]).reset_index(drop=True)
	acc_grouped_df["Experiment_target_lab"] = acc_grouped_df["Experiment target"].astype(str) + "_" + acc_grouped_df["Lab"].apply(lambda x : x.split(" ")[1].replace(",", ""))

	### Add Control accession for each control fastq accession using Rest API query:
	control_list = []
	for fastq_acc in acc_grouped_df.loc[:,"Controlled by"].tolist():
		URL = "https://www.encodeproject.org/biosample/" + fastq_acc + "/?frame=embedded" # frame=object; for short, brief and robust description
		response = requests.get(URL, headers=HEADERS) # GET the object
		response_json_dict = response.json() # Print the object; print json.dumps(response_json_dict, indent=4, separators=(',', ': ')) 		
		try:
			expt_control_acc = response_json_dict["dataset"].split("/")[-2] # extracts the experiment accession from the fastq accession:       
			control_list.append(expt_control_acc)
			#print expt_control_acc
		except KeyError: 
			print "Data not availabe to public for %s \n\n\n" %(fastq_acc_new)

	### For check with Encode portal - if needed:
	df_combined = pd.concat([pd.Series(control_list), acc_grouped_df.loc[:,"Controlled by"]], axis=1)

	### Give unique accession ID to each experiment and control:
	acc_grouped_df["Control"] = control_list
	acc_grouped_df["Control_file_acc_Id"] = acc_grouped_df["Control"].astype(str) + "_" + acc_grouped_df["Controlled by"].apply(lambda x : x.split(" ")[0])
	acc_grouped_df["Experiment_file_acc_Id"] = acc_grouped_df[["Experiment_acc_repId", "File accession"]].astype(str).apply(lambda x: '_'.join(x), axis=1)

	### Split into 2 reps:
	rep1_df = acc_grouped_df.iloc[::2]
	rep2_df = acc_grouped_df.iloc[1::2]

	if rep1_df.shape == rep2_df.shape:
		print("\n\nRep1 and Rep2 dataframe are of equal sizes ....\n")
	else:
		print("\nRecheck the rep1 and rep2 dataframe, aren't of equal sizes....\n")

	"""Create experiment accession libraries list for peak calling step"""
	# Rep1 and ctrl1 | Rep2 and ctrl2
	rep1_list_df = rep1_df.loc[:,"Experiment_acc_repId"]
	control1_list_df = rep1_df.loc[:,"Control_file_acc_Id"]
	rep2_list_df = rep2_df.loc[:,"Experiment_acc_repId"]
	control2_list_df = rep2_df.loc[:,"Control_file_acc_Id"]

	### Save sub-files for ENCODE pipeline processing:
	rep1_list_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_tf_rep1.txt"), sep="\t", header=None, index=False)
	rep2_list_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_tf_rep2.txt"), sep="\t", header=None, index=False)
	control1_list_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_tf_control1.txt"), sep="\t", header=None, index=False)
	control2_list_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_tf_control2.txt"), sep="\t", header=None, index=False)

	### Combine replicates-control info for later reference:
	rep1_list_df_reset = rep1_list_df.reset_index(drop=True)
	rep2_list_df_reset = rep2_list_df.reset_index(drop=True)
	control1_list_df_reset = control1_list_df.reset_index(drop=True)
	control2_list_df_reset = control2_list_df.reset_index(drop=True)
	combined_df = pd.concat([rep1_list_df_reset, rep2_list_df_reset, control1_list_df_reset, control2_list_df_reset], axis=1)
	combined_df.columns = ["Experiment_acc_repId", "Experiment_acc_repId_2", "Control_file_acc_Id", "Control_file_acc_Id_2"]
	combined_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_replicates_control_reference_file.txt"), sep="\t", header=True, index=True)

	### Create symbolic link or check if file was downloaded :
	downloaded_file_list = [ basename(each_fastq).split(".")[0] for each_fastq in glob(join(fastq_source_lib_dir,"*fastq.gz")) ]
	dframe_fastq_list = pd.concat([acc_grouped_df["File accession"], acc_grouped_df["Controlled by"]], ignore_index =True)
	to_download_list = [ fastq_id for fastq_id in dframe_fastq_list if fastq_id not in downloaded_file_list ]
	if len(to_download_list) > 0:
		print "\nReminder: (Fastq that needs to be downloaded)\n", to_download_list
		#os.system("wget --directory-prefix=fastq_dest_lib_dir https://www.encodeproject.org/files/${fastq_id}/@@download/${fastq_id}.fastq.gz")

	acc_dframe_colname_list = ["experiment_library", "experiment_fastq", "control_library", "control_fastq", "lab"]
	acc_dframe_value_list = []

	### Iterate via dframe index to structure each expts/controls lib dir and fastqs:
	acc_index_df = acc_grouped_df.reset_index(drop=True)
	for index in acc_index_df.index.values:
		experiment_library = acc_index_df["Experiment_acc_repId"].iloc[index]
		experiment_fastq = acc_index_df["File accession"].iloc[index] + ".fastq.gz"
		control_library = acc_index_df["Control_file_acc_Id"].iloc[index]
		control_fastq = acc_index_df["Controlled by"].iloc[index] + ".fastq.gz"
		lab = acc_index_df["Experiment_target_lab"].iloc[index] 
		acc_dframe_value_list.append([experiment_library, experiment_fastq, control_library, control_fastq, lab])
		#print index, experiment_library, experiment_fastq, control_library, control_fastq, lab
		
		dest_expt_lib_dir = join(fastq_dest_lib_dir, experiment_library)
		dest_control_lib_dir = join(fastq_dest_lib_dir, control_library)	
		
		source_expt_file = join(fastq_source_lib_dir, experiment_fastq)
		source_control_file = join(fastq_source_lib_dir, control_fastq)
		dest_expt_file = join(dest_expt_lib_dir, experiment_fastq)
		dest_control_file = join(dest_control_lib_dir, control_fastq)

		# if not exists(dest_expt_lib_dir):
		# 	makedirs(dest_expt_lib_dir)

		# if not exists(dest_control_lib_dir):
		# 	makedirs(dest_control_lib_dir) 

		# if not exists(dest_expt_file):
		# 	os.system("ln -fs %s %s" %(source_expt_file, dest_expt_lib_dir))
			
		# if not exists(dest_control_file):
		# 	os.system("ln -fs %s %s" %(source_control_file, dest_control_lib_dir))

	acc_dframe_dict = dict(zip(acc_dframe_colname_list, zip(*acc_dframe_value_list)))
	file_reference_dframe = pd.DataFrame(acc_dframe_dict)
	file_reference_dframe.to_csv(join(se_liblist_outdir, cell_type+"_SE_fastq_reference_metadata_final.txt"), sep="\t", header=True, index=True)
	return(file_reference_dframe, to_download_list)

se_reference_metadata, se_to_download = structure_SE_fastq_libraries_for_peakcall(metadata_tsv_file, cell_type, fastq_source_lib_dir, fastq_dest_lib_dir)



#cell_type = "GM12878"
#control_metadata_tsv_file = "control_Blood_chipseq_metadata_final.tsv"
cell_type = "K562"
metadata_tsv_file = "Blood_chipseq_metadata_final.tsv"
fastq_source_lib_dir = join(input_dir, "fastqs")
fastq_dest_lib_dir = join(input_dir, "k562_PE_fastq_libraries")

def structure_PE_fastq_libraries_for_peakcall(metadata_tsv_file, cell_type, fastq_source_lib_dir, fastq_dest_lib_dir):
	if not exists(fastq_dest_lib_dir):
		makedirs(fastq_dest_lib_dir)

	df = pd.read_csv(join(input_dir, metadata_tsv_file), sep="\t")
	""" Filter Ended fastq files """
	#df_se = df[df.loc[:, "Run type"] == "single-ended"]
	df_pe = df[df.loc[:, "Run type"] == "paired-ended"]
	final_df = df_pe[(df_pe.loc[:, "File format"] == "fastq") & (df_pe.loc[:, "File Status"] != "archived")]
	select_cols = ["File accession", 'Controlled by', "Experiment accession", "Biological replicate(s)", "Lab",
						"Paired with", "Paired end", "Experiment target", "Biosample term name"]
	selected_final_df = final_df.loc[:, select_cols]
	celltype_df = selected_final_df[selected_final_df.loc[:,"Biosample term name"]== cell_type]
	celltype_df_wo_eGFP = celltype_df[~celltype_df.loc[:,"Experiment target"].str.contains("eGFP")]
	
	select_cols = ['File accession', "Experiment target", 'Controlled by', 'Paired end','Experiment accession', 'Biological replicate(s)', 'Paired with']
	sorted_dframe = celltype_df_wo_eGFP.sort_values(["Experiment accession", "Biological replicate(s)"])
	sorted_dframe = sorted_dframe.reset_index(drop=True)

	### Select the samples with only controls available:
	sorted_dframe = sorted_dframe[~sorted_dframe['Controlled by'].isnull()]
	sorted_dframe["Controlled by"] = sorted_dframe['Controlled by'].replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True) # capture regex

	### Filter out sample with 1 reps or selecting samples with paired reps only:
	rep_idx = []
	expt_acc_list = sorted_dframe.loc[:,"Experiment accession"].tolist()
	for i,item in enumerate(expt_acc_list):
		if (expt_acc_list.count(item) >= 4):
			#print(item), i 
			rep_idx.append(i)

	final_rep_df = sorted_dframe.iloc[rep_idx]
	final_rep_df = final_rep_df.reset_index(drop=True)

	### Select top(head-2) replicates from the list (if more than 2 reps for any experiment):
	acc_grouped_df = final_rep_df.groupby(["Experiment accession"]).head(4).reset_index(drop=True)
	acc_grouped_df["Experiment_acc_repId"] = final_rep_df[["Experiment accession", "Biological replicate(s)"]].astype(str).apply(lambda x: '_'.join(x), axis=1)
	acc_grouped_df = acc_grouped_df.sort_values(["Experiment_acc_repId"]).reset_index(drop=True)
	acc_grouped_df["Experiment_target_lab"] = acc_grouped_df["Experiment target"].astype(str) + "_" + acc_grouped_df["Lab"].apply(lambda x : x.split(" ")[1].replace(",", ""))
	acc_grouped_df["Paired end"] = acc_grouped_df.loc[:,"Paired end"].map(int)
	acc_grouped_df["File_acc_pairId"] = acc_grouped_df[["File accession", "Paired end"]].astype(str).apply(lambda x: '_'.join(x), axis=1)
	acc_grouped_df["File_acc_control_pairId"] = acc_grouped_df[["Controlled by", "Paired end"]].astype(str).apply(lambda x: '_'.join(x), axis=1)
	
	### Add Control accession for each control fastq accession using Rest API query:
	control_list = []
	for fastq_acc in acc_grouped_df.loc[:,"Controlled by"].tolist():
		URL = "https://www.encodeproject.org/biosample/" + fastq_acc + "/?frame=embedded" # frame=object; for short, brief and robust description
		response = requests.get(URL, headers=HEADERS) # GET the object
		response_json_dict = response.json() # Print the object; print json.dumps(response_json_dict, indent=4, separators=(',', ': ')) 		
		try:
			#expt_control_acc_ENCFF = response_json_dict["dataset"].split("/")[-2] # extracts the experiment accession from the fastq accession:       
			expt_control_acc_ENCLB = response_json_dict[u'replicate'][u'library'][u'accession'] # extracts the experiment accession from the fastq accession:       
			control_list.append(expt_control_acc_ENCLB)
			print expt_control_acc_ENCLB
		except KeyError: 
			print "Data not availabe to public for %s \n\n\n" %(fastq_acc_new)

	### For check with Encode portal - if needed:
	df_combined = pd.concat([pd.Series(control_list), acc_grouped_df.loc[:,"Controlled by"]], axis=1)

	### Give unique accession ID to each experiment and control:
	acc_grouped_df["Control_ENCLB_accession"] = control_list
	acc_grouped_df["Control_file_acc_Id"] = acc_grouped_df["Control_ENCLB_accession"].astype(str) + "_" + acc_grouped_df["Controlled by"].apply(lambda x : x.split(" ")[0])
	acc_grouped_df["Experiment_file_acc_Id"] = acc_grouped_df[["Experiment_acc_repId", "File accession"]].astype(str).apply(lambda x: '_'.join(x), axis=1)

	### Split into 2 reps:
	acc_grouped_df = acc_grouped_df.sort_values("Experiment_acc_repId").reset_index(drop=True)
	rep1_df = acc_grouped_df.groupby(["Experiment accession"]).head(2).reset_index(drop=True)
	rep2_df = acc_grouped_df.groupby(["Experiment accession"]).tail(2).reset_index(drop=True)

	if rep1_df.shape == rep2_df.shape:
		print("\n\nRep1 and Rep2 dataframe are of equal sizes ....\n")
	else:
		print("\nRecheck the rep1 and rep2 dataframe, aren't of equal sizes....\n")

	"""Create experiment accession libraries list for peak calling step"""
	# Rep1 and ctrl1 | Rep2 and ctrl2
	rep1_list_df = rep1_df.loc[:,"Experiment_acc_repId"]
	control1_list_df = rep1_df.loc[:,"Control_ENCLB_accession"]
	rep2_list_df = rep2_df.loc[:,"Experiment_acc_repId"]
	control2_list_df = rep2_df.loc[:,"Control_ENCLB_accession"]

	""" Filter out the duplicates due to PE data nature to avoid peak calling of same replicate-control combination"""
	#combined_df = pd.concat([pd.Series(rep1_list_df), pd.Series(rep2_list_df), pd.Series(control1_list_df), pd.Series(control2_list_df)], axis=1)
	combined_df = pd.concat([rep1_list_df, rep2_list_df, control1_list_df, control2_list_df], axis=1)
	combined_df.columns = ["Experiment_acc_repId", "Experiment_acc_repId_2", "Control_ENCLB_accession", "Control_ENCLB_accession_2"]
	combined_df["new"] = combined_df.apply(lambda x: ":".join(x), axis=1)
	combined_df_list = combined_df["new"].tolist()
	uniq_df_list = list(set(combined_df_list))
	uniq_dframe_list = [ each.split(":") for each in uniq_df_list ]

	header=["Experiment_acc_repId", "Experiment_acc_repId_2", "Control_ENCLB_accession", "Control_ENCLB_accession_2"]
	filtered_uniq_df = pd.DataFrame(dict(zip(header, zip(*uniq_dframe_list))))
	filtered_uniq_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_replicates_control_reference_file.txt"), sep="\t", header=True, index=True)

	rep1_list_df = filtered_uniq_df.loc[:,"Experiment_acc_repId"]
	rep2_list_df = filtered_uniq_df.loc[:,"Experiment_acc_repId_2"]
	control1_list_df = filtered_uniq_df.loc[:,"Control_ENCLB_accession"]
	control2_list_df = filtered_uniq_df.loc[:,"Control_ENCLB_accession_2"]

	### Save sub-files for ENCODE pipeline processing:
	rep1_list_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_tf_rep1.txt"), sep="\t", header=None, index=False)
	rep2_list_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_tf_rep2.txt"), sep="\t", header=None, index=False)
	control1_list_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_tf_control1.txt"), sep="\t", header=None, index=False)
	control2_list_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_tf_control2.txt"), sep="\t", header=None, index=False)

	### Create symbolic link or check if file was downloaded :
	downloaded_file_list = [ basename(each_fastq).split(".")[0] for each_fastq in glob(join(fastq_source_lib_dir,"*fastq.gz")) ]
	dframe_fastq_list = pd.concat([acc_grouped_df["File accession"], acc_grouped_df["Controlled by"]], ignore_index =True)
	to_download_list = [ fastq_id for fastq_id in dframe_fastq_list if fastq_id not in downloaded_file_list ]
	if len(to_download_list) > 0:
		print "\nReminder: (Fastq that needs to be downloaded)\n", to_download_list
		#os.system("wget --directory-prefix=fastq_dest_lib_dir https://www.encodeproject.org/files/${fastq_id}/@@download/${fastq_id}.fastq.gz")

	acc_dframe_colname_list = ["experiment_library", "experiment_fastq", "paired_mate", "control_library", "control_fastq", "lab"]
	acc_dframe_value_list = []

	### Iterate via dframe index to structure each expts/controls lib dir and fastqs:
	acc_index_df = acc_grouped_df.reset_index(drop=True)
	for index in acc_index_df.index.values:
		experiment_library = acc_index_df["Experiment_acc_repId"].iloc[index]
		experiment_fastq = acc_index_df["File accession"].iloc[index] + ".fastq.gz"
		control_library = acc_index_df["Control_ENCLB_accession"].iloc[index]
		control_fastq = acc_index_df["Controlled by"].iloc[index] + ".fastq.gz"
		lab = acc_index_df["Experiment_target_lab"].iloc[index]
		paired_mate = acc_index_df["Paired end"].iloc[index].astype(str)

		new_experiment_fastq = acc_index_df["File accession"].iloc[index] + "_" + paired_mate + ".fastq.gz"
		new_control_fastq = acc_index_df["Controlled by"].iloc[index] + "_" + paired_mate + ".fastq.gz"
		acc_dframe_value_list.append([experiment_library, experiment_fastq, paired_mate, control_library, control_fastq, lab])
		#print index, experiment_library, experiment_fastq, control_library, control_fastq, lab
		
		dest_expt_lib_dir = join(fastq_dest_lib_dir, experiment_library)
		dest_control_lib_dir = join(fastq_dest_lib_dir, control_library)	
		
		source_expt_file = join(fastq_source_lib_dir, experiment_fastq)
		source_control_file = join(fastq_source_lib_dir, control_fastq)
		dest_expt_file = join(dest_expt_lib_dir, new_experiment_fastq)
		dest_control_file = join(dest_control_lib_dir, new_control_fastq)

		# if not exists(dest_expt_lib_dir):
		# 	makedirs(dest_expt_lib_dir)

		# if not exists(dest_control_lib_dir):
		# 	makedirs(dest_control_lib_dir) 

		# if not exists(dest_expt_file):
		# 	os.system("ln -fs %s %s" %(source_expt_file, dest_expt_file))
			
		# if not exists(dest_control_file):
		# 	os.system("ln -fs %s %s" %(source_control_file, dest_control_file))

	acc_dframe_dict = dict(zip(acc_dframe_colname_list, zip(*acc_dframe_value_list)))
	file_reference_dframe = pd.DataFrame(acc_dframe_dict)
	file_reference_dframe.to_csv(join(pe_liblist_outdir, cell_type+"_PE_fastq_reference_metadata_final.txt"), sep="\t", header=True, index=True)
	return(file_reference_dframe, to_download_list)

pe_reference_metadata, pe_to_download = structure_PE_fastq_libraries_for_peakcall(metadata_tsv_file, cell_type, fastq_source_lib_dir, fastq_dest_lib_dir)


### For single-ended files:
input_file1 = join(se_liblist_outdir, "K562_SE_replicates_control_reference_file.txt")
input_file2 = join(se_liblist_outdir, "K562_SE_fastq_reference_metadata_final.txt")
run_type = "SE"
cell_type = "K562"

def merge_SE_files_accession_for_labinfo(input_file1, input_file2, run_type, cell_type):
	df1 = pd.read_csv(input_file1, sep = "\t")
	df1_select = df1.iloc[:, [1, 2, 3, 4]]
	df1_select["Run_type"] = run_type 
	df1_select.columns = ["Expt_rep", "Expt_rep2", "Control_rep", "Control_rep2", "Run_type"]

	df2 = pd.read_csv(input_file2, sep = "\t")
	df2_select = df2.iloc[:, [4,5]]
	df2_select.columns = ["Expt_rep", "lab"]

	all_dfs = [df1_select, df2_select]
	merged_df = reduce(lambda left, right: pd.merge(left, right, on=["Expt_rep"]), all_dfs)
	merged_df.to_csv(join(se_liblist_outdir, cell_type+"_SE_replicates_control_reference_file_with_labinfo.txt"), sep="\t", header=True, index=True)
	return(merged_df)

se_merged_df = merge_SE_files_accession_for_labinfo(input_file1, input_file2, run_type, cell_type)

### For pair-ended files:
input_file1 = join(pe_liblist_outdir, "K562_PE_replicates_control_reference_file.txt")
input_file2 = join(pe_liblist_outdir, "K562_PE_fastq_reference_metadata_final.txt")
run_type = "PE"
cell_type = "K562"

def merge_PE_files_accession_for_labinfo(input_file1, input_file2, run_type, cell_type):
	df1 = pd.read_csv(input_file1, sep = "\t")
	df1_select = df1.iloc[:, [3, 4, 1, 2]]
	df1_select["Run_type"] = run_type 
	df1_select.columns = ["Expt_rep", "Expt_rep2", "Control_rep", "Control_rep2", "Run_type"]

	df2 = pd.read_csv(input_file2, sep = "\t")
	df2_select = df2.iloc[:, [4,5]]
	df2_select = df2_select.drop_duplicates(["experiment_library"])
	df2_select.columns = ["Expt_rep", "lab"]

	all_dfs = [df1_select, df2_select]
	merged_df = reduce(lambda left, right: pd.merge(left, right, on=["Expt_rep"]), all_dfs)
	merged_df.to_csv(join(pe_liblist_outdir, cell_type+"_PE_replicates_control_reference_file_with_labinfo.txt"), sep="\t", header=True, index=True)
	return(merged_df)

pe_merged_df = merge_PE_files_accession_for_labinfo(input_file1, input_file2, run_type, cell_type)

final_combined_df = pd.concat([se_merged_df, pe_merged_df], ignore_index=True)
final_combined_df.to_csv(join(output_dir, cell_type + "_replicates_control_reference_file_with_labinfo_combined.txt"), sep="\t", header=True, index=True)


tf_file1_k562 = os.path.expanduser("~/Dropbox/for_encode/encode3_paper/chipseq_metadata_info/K562_replicates_control_reference_file_with_labinfo_combined.txt")
tf_file2_hepg2 = os.path.expanduser("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx")
output_dir = os.path.expanduser("~/Dropbox/for_encode/encode3_paper/chipseq_metadata_info")

def find_tf_overlap_subsets(tf_file1, tf_file2):
	tf_file1 = tf_file1_k562
	df_k562 = pd.read_csv(tf_file1, sep="\t")
	df_k562["TF_name"] = df_k562["lab"].apply(lambda x : x.replace("FLAG-", "").split("-")[0]).str.upper()	
	k562_tf_list = df_k562["TF_name"].sort_values().unique()

	tf_file2 = tf_file2_hepg2
	xls = pd.ExcelFile(excel_file)
	df_hepg2 = xls.parse("unique_TFs")
	df_hepg2["TF_name"] = df_hepg2["Target"].apply(lambda x : x.replace("[FLAG]", "").split("_")[0]).str.upper()
	hepg2_tf_list = df_hepg2["TF_name"].sort_values().unique()
	
	df1 = df_hepg2["TF_name"].sort_values().drop_duplicates().reset_index(drop=True)
	df2 = df_k562["TF_name"].sort_values().drop_duplicates().reset_index(drop=True)	
	df_overlap = pd.merge(pd.DataFrame(df1), pd.DataFrame(df2), on="TF_name")
	#pd.concat([df_overlap1, df_overlap2], axis=1)

	hepg2_overlap_list = hepg2_tf_list[pd.Series(hepg2_tf_list).isin(k562_tf_list)]
	hepg2_non_overlap_list = hepg2_tf_list[~pd.Series(hepg2_tf_list).isin(k562_tf_list)]
	pd.Series(hepg2_overlap_list).to_csv(join(output_dir, "HepG2_overlapping_tfs_file1.txt"), sep="\t", header=True, index=True)
	pd.Series(hepg2_non_overlap_list).to_csv(join(output_dir, "HepG2_non_overlapping_tfs_file1.txt"), sep="\t", header=True, index=True)

	k562_overlap_list = k562_tf_list[pd.Series(k562_tf_list).isin(hepg2_tf_list)]
	k562_non_overlap_list = k562_tf_list[~pd.Series(k562_tf_list).isin(hepg2_tf_list)]
	pd.Series(k562_non_overlap_list).to_csv(join(output_dir, "K562_non_overlapping_tfs_file2.txt"), sep="\t", header=True, index=True)
	return(hepg2_overlap_list)

overlapped_tfs = find_tf_overlap_subsets(tf_file1_k562, tf_file2_hepg2)



### Step 2: Find the number of peaks from all experiments and select the one with highest number of peaks for
### overlapping TFs or experiments.

#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import shutil
from glob import glob
from os.path import join
from os.path import basename
from os.path import expanduser

reference_metadata_file = expanduser("~/k562_chip/chipseq_metadata_info/K562_replicates_control_reference_file_with_labinfo_combined.txt")
metadata_file_dir = expanduser("~/k562_chip/chipseq_metadata_info")

idr_peaks_dir = expanduser("~/k562_chip/idr_passed_peaks_total")
uniq_tf_dir = join(idr_peaks_dir, "unique_TFs")
if not os.path.exists(uniq_tf_dir):
    os.makedirs(uniq_tf_dir)

def compute_peak_count_with_unique_TFs(reference_metadata_file, idr_peaks_dir, metadata_file_dir ):
	df = pd.read_csv(reference_metadata_file, sep="\t")
	peak_count_list = []
	tf_name_list = []
	lab_name_list = []
	for index, row in df.iterrows():
		SL_rep1 = row[u'Expt_rep']
		SL_rep2 = row[u'Expt_rep2']
		tf_name = row["lab"].replace("FLAG-","").split("-")[0]
		lab_name = row[u'lab'].replace("FLAG-","").split("-")[1].split("_")[1]
		tf_name_list.append(tf_name)
		lab_name_list.append(lab_name)
		print index, SL_rep1, SL_rep2, tf_name, lab_name

		try:
			full_path = expanduser(join(idr_peaks_dir, SL_rep1 + "*" + SL_rep2 + "*"))
			file_name = glob(full_path)[0]
			#print file_name
			num_lines = sum(1 for line in open(file_name))
			print tf_name, num_lines, lab_name, "\n"
			peak_count_list.append(num_lines)

		except IndexError:
			#print file_name
			num_lines = 0
			print tf_name, num_lines, lab_name, "\n"
			peak_count_list.append(num_lines)

	### maintains the same order as the original df file:
	df["Peak_count"] = peak_count_list
	df["TF_name"] = tf_name_list
	df["Lab_name"] = lab_name_list
	df_grouped = df.iloc[df.groupby(["TF_name"])["Peak_count"].idxmax()]
	df_sorted = df_grouped.sort_values(["TF_name"])
	
	### Files with highest peak count:
	df_final1 = df_sorted.loc[:,["Expt_rep", "Expt_rep2", "TF_name", "Lab_name", "Run_type", "Peak_count"]]
	df_final2 = df_sorted.loc[:,["TF_name", "Lab_name", "Run_type", "Peak_count"]]
	df_final1.to_csv(join(metadata_file_dir, "Encode_k562_full_datasets_with_peakcount"), header=True, sep="\t", index=False)
	df_final1.to_excel(join(metadata_file_dir, "Encode_k562_full_datasets_with_peakcount.xls"), header=True, sheet_name="Sheet1")
	df_final2.to_csv(join(metadata_file_dir, "Encode_k562_datasets_with_peakcount"), header=True, sep="\t", index=False)
	df_final2.to_excel(join(metadata_file_dir, "Encode_k562_datasets_with_peakcount.xls"), header=True, sheet_name="Sheet1")
	return(df_final2)

unique_tfs_peak_count = compute_peak_count_with_unique_TFs(reference_metadata_file, idr_peaks_dir, output_dir )


metadata_copyfrom_file = join(metadata_file_dir, "Encode_k562_full_datasets_with_peakcount")
def copy_unique_TF_files_to_dir(metadata_copyfrom_file, idr_peaks_dir, uniq_tf_dir):
	df = pd.read_csv(metadata_copyfrom_file, sep="\t")
	df_final = df.loc[df["Peak_count"] > 0 ]
	for index, row in df_final.iterrows():
		SL_rep1 = row[u'Expt_rep']
		SL_rep2 = row[u'Expt_rep2']
		tf_name = row["TF_name"]
		full_path = os.path.expanduser(join(idr_peaks_dir, SL_rep1 + "*" + SL_rep2 + "*" ))
		source_file_name = glob(full_path)[0] #print file_name
		dest_file_name = "ENCODE_K562.nodup.IDR0.02.filt.narrowPeak_" + tf_name  #print file_name
		#dest_file_name = basename(glob(full_path)[0]) + "_" + tf_name #print file_name
		shutil.copy(source_file_name, join(uniq_tf_dir,dest_file_name))
		#os.system("cp %s %s" %(source_file_name, join(uniq_tf_dir,dest_file_name)))
	return()

file_transfer = copy_unique_TF_files_to_dir(metadata_copyfrom_file, idr_peaks_dir, uniq_tf_dir)


from reportlab.pdfgen import canvas
from os.path import basename, splitext, join

output_dir = "/Users/suryachhetri/Dropbox/for_encode/encode3_paper/chipseq_metadata_info/hepg2_tf"
file_list = [each for each in glob("*.pdf")]
for each in file_list:
    filename = splitext(basename(each))[0]
    c = canvas.Canvas(join(output_dir, filename + "_metadata" + ".pdf"))
    c.drawString(100,750,"Welcome to " + filename)
    c.save()



