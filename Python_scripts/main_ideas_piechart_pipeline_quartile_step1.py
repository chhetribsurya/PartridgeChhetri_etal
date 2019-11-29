#!/usr/bin/env python

import glob
import os
import re
from os.path import join
from os.path import basename

bsub_options = "bsub -We 10:00 -q c7normal -n 1 -R span[hosts=1]"

# Input dir path containing files:
dir_path = "/gpfs/gpfs1/home/schhetri/paper_revision_3/quartile_motif_analysis/tf_files/quartile_tf_files/*quartile*.bed"
file_list = glob.glob(dir_path)

# Input target output dir:
output_dir = os.path.expanduser("~/paper_revision_3/piecharts_quartile_motifs")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#regex = re.compile(r".*narrowPeak_(.*)$")
regex = re.compile(r".*/(.*).srt.(.*).bed")
for each_tf_file in file_list:
    capture_list = regex.findall(each_tf_file)[0]
    tf_name = capture_list[0] + "." + capture_list[1]
    # general environment variables, preparing for pipeline run:   
    os.environ["peak_file_full_path"] = each_tf_file
    os.environ["bsub_options"] = bsub_options
    os.environ["log_name"] = "quartile_tf_overlap_piechart.out"
    os.environ["log_dir"] = "."
    os.environ["bjob_name"] = "Intersection and ideas piechart for " + tf_name.split("[")[0] # Because of flag suffix with TFs; HLF[FLAG].split("[")[0]
    os.environ["output_dir"] = output_dir 
    os.environ["piechart_calling_script"] = "/gpfs/gpfs1/home/schhetri/paper_revision_3/ideas_piechart_hepg2_total_quartile.py" 
    os.system('$bsub_options -J "$bjob_name" -o $log_dir/$log_name "python $piechart_calling_script $peak_file_full_path $output_dir"')

