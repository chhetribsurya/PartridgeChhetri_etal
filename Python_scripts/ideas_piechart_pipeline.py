#!/usr/bin/env python

import glob
import os
import re
from os.path import join
from os.path import basename

bsub_options = "bsub -We 10:00 -q normal -n 1 -R span[hosts=1]"

# Input dir path containing files:
dir_path = os.path.expanduser("~/batch_I/idr_passed_peaks_total/unique_TFs/SL*") 
file_list = glob.glob(dir_path)

# Input target output dir:
output_dir = os.path.expanduser("~/batch_I/piecharts_ideas_unique")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

regex = re.compile(r".*narrowPeak_(.*)$")
for each_tf_file in file_list:
    tf_name = regex.findall(basename(each_tf_file))[0]
    # general environment variables, preparing for pipeline run:   
    os.environ["peak_file_full_path"] = each_tf_file
    os.environ["bsub_options"] = bsub_options
    os.environ["log_name"] = "tf_overlap_piechart.out"
    os.environ["log_dir"] = "."
    os.environ["bjob_name"] = "Intersection and ideas piechart for " + tf_name.split("[")[0] # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["output_dir"] = output_dir 
    os.environ["piechart_calling_script"] = "./batch_I/ideas_piechart_hepg2_total.py" 
    os.system('$bsub_options -J "$bjob_name" -o $log_dir/$log_name "python $piechart_calling_script $peak_file_full_path $output_dir"')

