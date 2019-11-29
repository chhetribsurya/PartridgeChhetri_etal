#!/usr/bin/env python
import os, sys, re
import pandas as pd, numpy as np
import pybedtools
import json, time, subprocess
from os.path import basename, splitext, join
import glob

start_time = time.time()
full_file_path = sys.argv[1]
output_dir = sys.argv[2]

file_suffix = "intersectbed_counts_with_ideas.txt"

####Input files for script to run are attached in the dir input_file_for script:
ideas_hepg2_file = os.path.expanduser("~/batch_I/hepg2_ideas_36_dense.bed")  

#Provide the file path location to save the output of ENCODE_TF_cluster file as JSON format (could be useful for future usage):
JSON_dict_file = os.path.expanduser(os.path.join(output_dir,"ENCODE_ideas_JSON_1.txt"))

def parse_ideas(ideas_file, Json_output_file):
	with open(ideas_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file.readlines():
				#print line
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, ch_state = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(ch_state)]

				if ch_state in master_dict_return.keys():
					if chrom in master_dict_return[ch_state].keys():
						master_dict_return[ch_state][chrom].append((chrom, int(start), int(end), "\t".join(line_info)))
					else:
						master_dict_return[ch_state][chrom] = [(chrom, int(start), int(end), "\t".join(line_info))]
				else:
					master_dict_return[ch_state] = {chrom:[(chrom, int(start), int(end), "\t".join(line_info))]}
					master_dict_return[ch_state].update({"A_bed_count_hits":0})
					master_dict_return[ch_state].update({"B_bed_count_hits":0})
					master_dict_return[ch_state].update({"custom_overlap_list":[]})
					master_dict_return[ch_state].update({"As_overlap_list":[]})
					master_dict_return[ch_state].update({"Bs_overlap_list":[]})

			json.dump(master_dict_return, outfile)
			#print master_dict_return
			return(master_dict_return)

#master_ideas_dict = parse_ideas(ideas_hepg2_file, JSON_dict_file)
#print "\n\nTime to parse ideas_segmentation based on region and chromosome = ", time.time()-start_time

def parse_ideas_based_on_region(ideas_file, Json_output_file):
	with open(ideas_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file.readlines():
				#print line
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, ch_state = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(ch_state)]

				if ch_state in master_dict_return.keys():
					master_dict_return[ch_state]["loci"].append([chrom, int(start), int(end)])
					#master_dict_return[ch_state]["loci"].append([chrom, int(start), int(end), "\t".join(line_info)])
					
				else:
					master_dict_return[ch_state] = {"loci" : [[chrom, int(start), int(end)]]}
					#master_dict_return[ch_state] = {"loci" : [[chrom, int(start), int(end), "\t".join(line_info)]]}
					master_dict_return[ch_state].update({"A_bed_count_hits":0})
					master_dict_return[ch_state].update({"B_bed_count_hits":0})
					master_dict_return[ch_state].update({"custom_overlap_list":[]})
					master_dict_return[ch_state].update({"As_overlap_list":[]})
					master_dict_return[ch_state].update({"Bs_overlap_list":[]})

			json.dump(master_dict_return, outfile)
			#print master_dict_return
			return(master_dict_return)

def load_stringf_for_peak_files(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df = df.iloc[:,0:6]
	df.columns = ["chrom", "start", "end", "name", "score", "strand"]
	TF_peak_file = df.loc[:,["chrom","start", "end"]]
	TF_peak_list =	TF_peak_file.values.tolist()
	TF_peak_string_list = [ "\t".join(map(str, line)) for line in TF_peak_list]
	TF_peak_string_bed_format = "\n".join(TF_peak_string_list)
	TF_bed_file = pybedtools.BedTool(TF_peak_string_bed_format,from_string=True)
	return(TF_bed_file)

def load_stringf_for_centerOfpeak_files(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df = df.iloc[:,0:6]
	df.columns = ["chrom", "start1", "end1", "name", "score", "strand"]
	df["mid_point"] = (df["start1"] + df["end1"])/2
	df["start"] = df["mid_point"].astype(int)
	df["end"] = df["mid_point"].astype(int) + 1

	TF_peak_file = df.loc[:,["chrom","start", "end"]]
	TF_peak_list =	TF_peak_file.values.tolist()
	TF_peak_string_list = [ "\t".join(map(str, line)) for line in TF_peak_list]
	TF_peak_string_bed_format = "\n".join(TF_peak_string_list)
	TF_bed_file = pybedtools.BedTool(TF_peak_string_bed_format,from_string=True)
	return(TF_bed_file)

tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(full_file_path))[0]
file_path = full_file_path
print "\nCurrently Processing %s\n" %(tf_name)
#final_output_file = output_dir + "/" + tf_name + "_" + "intersectbed_counts_with_ideas.txt"
final_output_file = os.path.join(output_dir, (tf_name + "_" + file_suffix))
if not os.path.exists(final_output_file):  
    peaks_bedfile1 = file_path
    master_dict_return = {}
    master_ideas_dict = parse_ideas_based_on_region(ideas_hepg2_file, JSON_dict_file)
    print "\n\nTime to parse ideas_segmentation based on region = ", time.time()-start_time
    start_time = time.time()
    print "\n\n Overlap analysis starting..."
    peak_file = load_stringf_for_centerOfpeak_files(peaks_bedfile1)

    for key,value in master_ideas_dict.iteritems():
        # name the output file with its keyword:
        file_name = key.replace("/","-") + ".txt"
        outfile_path = os.path.join(output_dir, file_name)
        each_tf_dict = master_ideas_dict[key]
        print "processing....",  key, "final......."
        TF_loci_list = each_tf_dict.get("loci")
        TF_loci_string_list = [ "\t".join(map(str, line)) for line in TF_loci_list]
        TF_loci_string_bed_format = "\n".join(TF_loci_string_list)
        ideas_loci_bed_file = pybedtools.BedTool(TF_loci_string_bed_format,from_string=True)
        bed_intersect = peak_file.intersect(ideas_loci_bed_file, wa = True, wb = True)
        intersect_count = bed_intersect.count()

        if intersect_count >= 1:
            # Bedtool object converted to pandas dataframe:
            df_bed = pd.read_table(bed_intersect.fn, sep="\t", header=None)
            As_overlap_list = df_bed.iloc[:,0:3]
            As_overlap_list.columns = ["chrom", "start", "end"]
            As_overlap_list = As_overlap_list.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
            Bs_overlap_list = df_bed.iloc[:,3:6]
            Bs_overlap_list.columns = ["chrom", "start", "end"]
            Bs_overlap_list = Bs_overlap_list.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
            custom_overlap_list = df_bed

            print key, "appending As overlap_list"
            master_ideas_dict[key]["As_overlap_list"].append(As_overlap_list)
            master_ideas_dict[key]["As_overlap_list"][0] = master_ideas_dict[key]["As_overlap_list"][0].drop_duplicates(["chrom", "start", "end"])
            master_ideas_dict[key]["A_bed_count_hits"] = len(master_ideas_dict[key]["As_overlap_list"][0].index)


            master_ideas_dict[key]["Bs_overlap_list"].append(Bs_overlap_list)
            master_ideas_dict[key]["Bs_overlap_list"][0] = master_ideas_dict[key]["Bs_overlap_list"][0].drop_duplicates(["chrom", "start", "end"])
            master_ideas_dict[key]["B_bed_count_hits"] = len(master_ideas_dict[key]["Bs_overlap_list"][0].index)

            master_ideas_dict[key]["custom_overlap_list"].append(df_bed)

            df_bed.to_csv(outfile_path, sep ="\t", header = True, index = False)
            outfile_path_prefix = os.path.splitext(outfile_path)
            master_ideas_dict[key]["As_overlap_list"][0].to_csv(outfile_path_prefix[0] + "_with_A_overlap.bed", sep="\t")

    print "\n overlap b/w bed files completed!!!!!"
    print "Time for overlap operation = ", time.time()-start_time

    out_file = output_dir + "/" + tf_name+ "_" + "intersectbed_counts_with_ideas.txt"
    with open(out_file, "w") as outfile:
        for key,value in master_ideas_dict.items():
            result = "%s\t%s" %(key, master_ideas_dict[key]["A_bed_count_hits"])
            outfile.write(result + "\n")
            print result

    print "\nIntersection job for %s completed....!!\n" %(tf_name)

#ideas_file = output_dir + "/" + tf_name+ "_" + "intersectbed_counts_with_ideas.txt"
ideas_file = os.path.join(output_dir, (tf_name + "_" + file_suffix))
os.environ["ideas_file"] =  ideas_file
os.environ["tf_name"] = tf_name
os.environ["out_dir_name"] = output_dir
os.system('Rscript ./batch_I/ideas_piechart_morgan.R $ideas_file $tf_name $out_dir_name')
#subprocess.call("Rscript " + "ideas_piechart_morgan.R " +  ideas_file + tf_name + output_dir], shell=True)
#subprocess.call("Rscript ideas_piechart_morgan.R --args ideas_file tf_name output_dir", shell=True)
print "\nRunning of Rscript for plotting completed!!!....\n"
print "\nCheck your plots in %s dir\n" %(output_dir)
	
	
    


    
    



    
    




