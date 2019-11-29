import os
import pandas as pd
from os.path import basename, join, splitext
from glob import glob

input_file = "/gpfs/gpfs1/home/schhetri/hek293/metadata/hek293_final_selected_metadata.txt"
file_dir = "/gpfs/gpfs1/home/schhetri/hek293/TF_files"
output_dir = "/gpfs/gpfs1/home/schhetri/hek293/TF_annotated_files"

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

### Download files in targeted directory (TF_files):
#for each_file in $(cat hek293_download_urls.txt); do wget $each_file -P ../TF_files/; done

### Create dictionary of ENCFF*-TF name and remove eGFP TF tags, for easy usage: 
df = pd.read_csv(input_file, sep="\t")
df["Experiment target"] = df["Experiment target"].apply(lambda x: x.replace("eGFP-", "").replace("-human", ""))
df_final = df.loc[:,["File accession", "Experiment target"]]
final_dict = df_final.set_index("File accession").T.to_dict("list")
files_list = glob(join(file_dir, "*.bed"))

for each_file in files_list:
	file_base = basename(each_file)
	encode_acc = splitext(file_base)[0]
	TF_name = final_dict[encode_acc][0] + "_" + encode_acc
	#os.environ["Output_dir"] = output_dir
	os.environ["File_name"] = each_file
	os.environ["TF_annotated_file"] = join(output_dir, TF_name + ".bed") 
	CMD = "cp $File_name $TF_annotated_file"
	os.system(CMD)

def combine_tf_bedfiles(TF_file_list_of_interest):
	tf_file_list = TF_file_list_of_interest
	concat_file_list = []
	for each_file in tf_file_list:
		tf_name = splitext(basename(each_file))[0]
		df_read = pd.read_csv(each_file, header=None, sep="\t")
		df_read = df_read.iloc[:,[0,1,2]]
		df_read.columns = ["peak_chrom", "start", "end"]
		df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
		df_read["peak_start"] = df_read["mid_point"].astype(int) - (50)
		df_read["peak_end"] = df_read["mid_point"].astype(int) + (50)
		df_read["strand"] = "."
		df_read["tf_name"] = tf_name
		select_cols = ["peak_chrom","peak_start","peak_end","start","end", "strand", "tf_name"]
		df_read = df_read.loc[:,select_cols]
		print "test\n\n"
		print df_read
		concat_file_list.append(df_read)

	combined_tf_df = pd.concat(concat_file_list, ignore_index=True)
	combined_tf_df.to_csv(join(output_dir, "final_TF_coord_df.bed"), header=False, index=False, sep="\t")
	return(combined_tf_df)

TF_annotated_list = glob(join(output_dir, "*.bed"))
combined_tf_df = combine_tf_bedfiles(TF_annotated_list) # cr_tf_file_list; dbf_tf_file_list; all_tf_file_list




