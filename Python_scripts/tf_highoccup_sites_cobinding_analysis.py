import pandas as pd
import pybedtools
import os
from os.path import join
from os.path import splitext
from os.path import basename


#main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_co_binding_total"
input_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_prom_enh_categorised/merged_d_100"
main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_hotspots_prom_enh_categorised/merged_d_100/cobinding_analysis"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

df_hotspot = pd.read_csv(join(input_dir, "final_tf_hotspots_sites_with_TF_details_100bp_test.bed"), sep="\t")
df_hotspot["count"] = 1
df_hotspot["tf_size"] = df_hotspot.groupby(["ideas_chrom", "ideas_start", "ideas_end", "state_anno"])["count"].transform(lambda x: x.sum())
df_50_125 = df_hotspot[(df_hotspot["tf_size"] >=50) & (df_hotspot["tf_size"] <=125)] 
df_grouped= df_50_125.groupby(["ideas_chrom", "ideas_start", "ideas_end", "state_anno"])

### Count the TFs occupancy on those TF hotpsot sites:
hotspot_tf_count = df_50_125.groupby(["tf_name"]).size().reset_index(name="count")
hotspot_tf_count.to_csv(join(main_dir, "tf_hotspots_cobinding_analysis_with_tf_counts.txt"), sep="\t", header=True, index=True)

### Break down the TFs occupancy on different regulatory region TF hotpsot sites:
hotspot_tf_reg_region_count = df_50_125.groupby(["tf_name", "state_anno"]).size().reset_index(name="count")
hotspot_tf_reg_region_count.to_csv(join(main_dir, "tf_hotspots_cobinding_analysis_with_region_based_tf_counts.txt"), sep="\t", header=True, index=True)

### Filter by rows of interest using .loc method:
whole_genome_hotspot = df_grouped.size().reset_index(name="count")
enh_assoc_hotspot = whole_genome_hotspot.loc[whole_genome_hotspot["state_anno"].isin(["sEnh_assoc", "wEnh_assoc"]),:]
prom_assoc_hotspot = whole_genome_hotspot.loc[whole_genome_hotspot["state_anno"]=="prom_assoc",:]


chrom_list = []
start_list = []
end_list = []
state_list = []

tf_cobind_list = []
tf_name_list = []

for group_name, data in df_grouped:
    data_tf_list = data["tf_name"].unique()
    for each_i, each in enumerate(data_tf_list):
        for other_i, other in enumerate(data_tf_list):
            if each_i != other_i:
                chrom_list.append(group_name[0])
				start_list.append(group_name[1])
				end_list.append(group_name[2])
				state_list.append(group_name[3])
				tf_name_list.append(each)
				tf_cobind_list.append(other)

new_hotspot_df = pd.DataFrame({"tf_name":tf_name_list, "cobind_tf_name" : tf_cobind_list, "chrom": chrom_list, "start": start_list, 
								"end" : end_list, "state_anno" : state_list})

### For whole genome tf_hotspots:
final_hotspot_df = new_hotspot_df.groupby(["tf_name", "cobind_tf_name"]).size().reset_index(name="size")
final_hotspot_df["total_sites"] = whole_genome_hotspot.shape[0]
final_hotspot_df["percent_overlap"] = (final_hotspot_df["size"]/final_hotspot_df["total_sites"])*100
final_hotspot_df.to_csv(join(main_dir, "tf_hotspots_cobinding_analysis_whole_genome.txt"), header=True, sep="\t")

final_hotspot_df_sorted = final_hotspot_df.sort_values(["percent_overlap"], ascending = False)
final_hotspot_df_sorted.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_whole_genome_sorted.txt"), sep="\t", header=True, index=True)

# Generating a heatmap data from the cobind peakcount dataframe:
final_hotspot_df_heatmap_table = final_hotspot_df_sorted.pivot(index="tf_name", columns="cobind_tf_name", values = "percent_overlap")
final_hotspot_df_heatmap_table.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_whole_genome_heatmap.txt"), sep="\t", header=True, index=True)

### For enhancer associated tf_hotspots:
enh_assoc_hotspot_df = new_hotspot_df.loc[new_hotspot_df["state_anno"].isin(["sEnh_assoc", "wEnh_assoc"]),:]
enh_assoc_hotspot_df = enh_assoc_hotspot_df.groupby(["tf_name", "cobind_tf_name"]).size().reset_index(name="size")
enh_assoc_hotspot_df["total_sites"] = enh_assoc_hotspot.shape[0]
enh_assoc_hotspot_df["percent_overlap"] = (enh_assoc_hotspot_df["size"]/enh_assoc_hotspot_df["total_sites"])*100
enh_assoc_hotspot_df.to_csv(join(main_dir, "tf_hotspots_cobinding_analysis_enh_assoc.txt"), header=True, sep="\t")

enh_assoc_hotspot_df_sorted = enh_assoc_hotspot_df.sort_values(["percent_overlap"], ascending = False)
enh_assoc_hotspot_df_sorted.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_enh_assoc_sorted.txt"), sep="\t", header=True, index=True)

# Generating a heatmap data from the cobind peakcount dataframe:
enh_assoc_hotspot_df_heatmap_table = enh_assoc_hotspot_df_sorted.pivot(index="tf_name", columns="cobind_tf_name", values = "percent_overlap")
enh_assoc_hotspot_df_heatmap_table.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_enh_assoc_heatmap.txt"), sep="\t", header=True, index=True)

### For promoter associated tf_hotspots:
prom_assoc_hotspot_df = new_hotspot_df.loc[new_hotspot_df["state_anno"]=="prom_assoc",:]
prom_assoc_hotspot_df = prom_assoc_hotspot_df.groupby(["tf_name", "cobind_tf_name"]).size().reset_index(name="size")
prom_assoc_hotspot_df["total_sites"] = prom_assoc_hotspot.shape[0]
prom_assoc_hotspot_df["percent_overlap"] = (prom_assoc_hotspot_df["size"]/prom_assoc_hotspot_df["total_sites"])*100
prom_assoc_hotspot_df.to_csv(join(main_dir, "tf_hotspots_cobinding_analysis_prom_assoc.txt"), header=True, sep="\t")

prom_assoc_hotspot_df_sorted = prom_assoc_hotspot_df.sort_values(["percent_overlap"], ascending = False)
prom_assoc_hotspot_df_sorted.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_prom_assoc_sorted.txt"), sep="\t", header=True, index=True)

# Generating a heatmap data from the cobind peakcount dataframe:
prom_assoc_hotspot_df_heatmap_table = prom_assoc_hotspot_df_sorted.pivot(index="tf_name", columns="cobind_tf_name", values = "percent_overlap")
prom_assoc_hotspot_df_heatmap_table.to_csv(join(main_dir,"tf_hotspots_cobinding_analysis_prom_assoc_heatmap.txt"), sep="\t", header=True, index=True)




