###########################
# Complex heatmap Library
###########################


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

output_dir <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/misc/TFs_Annotation_file.txt", sep="\t")
#anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))

# 2kb up and downstream heatmap:
input_file <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip/quartilebased_ideastyle_promoter_state_dist_for_heatmap.txt"

# input_file <- "/home/surya/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap_zscore.txt"
# input_file <- "~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt"

read_file <- fread(input_file, sep="\t", header=TRUE) 
read_file$Total_Frac <- NULL #delete last column
read_file$REST_STATE <- NULL
#read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged[,c("tf_name", "Category")]
#read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file)
data <- read_df[,2:ncol(read_df)] # Subtract last 4 columns
mat_data <- as.matrix(data)

# Rowwise scaling :As interested in binding variance across the cis-regulatory region
# for any given factor
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)

output_file_name <- file.path(output_dir, "Quartilebased_IDEAD_state_heatmap_clustering.pdf")       
#output_file_name <- file.path(output_dir, "Promoter_state_heatmap_figure_minus_rest.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS States",
    row_title="Quartile based Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 3.5),
    column_names_gp = gpar(fontsize = 6)
    # km=4,
    ## For heatmap and annotation legends (add new-legends)
    # heatmap_legend_param = list(
    # title_gp = gpar(fontsize = 10), 
    # labels_gp = gpar(fontsize = 6),
    # legend_direction = "horizontal",
    # legend_width = unit(3, "cm"), title_position = "topcenter"
    # )
) 

dev.off()

########################################
########################################

#Quartile plus 208 TFs based clustering 

########################################
########################################

output_dir <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/misc/TFs_Annotation_file.txt", sep="\t")

# 2kb up and downstream heatmap:
input_file <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip/quartilebased_plus208_ideastyle_promoter_state_dist_for_heatmap.txt"

read_file <- fread(input_file, sep="\t", header=TRUE) 
read_file$Total_Frac <- NULL #delete last column
read_file$REST_STATE <- NULL

read_df <- as.data.frame(read_file)
data <- read_df[,2:ncol(read_df)] # Subtract last 4 columns
mat_data <- as.matrix(data)

# Rowwise scaling :As interested in binding variance across the cis-regulatory region
# for any given factor
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)

output_file_name <- file.path(output_dir, "Quartilebased_plus208_IDEAS_state_heatmap_clustering.pdf")       
#output_file_name <- file.path(output_dir, "Promoter_state_heatmap_figure_minus_rest.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS States",
    row_title="Quartile based Factors plus 208",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.2),
    column_names_gp = gpar(fontsize = 6)
    # km=4,
    ## For heatmap and annotation legends (add new-legends)
    # heatmap_legend_param = list(
    # title_gp = gpar(fontsize = 10), 
    # labels_gp = gpar(fontsize = 6),
    # legend_direction = "horizontal",
    # legend_width = unit(3, "cm"), title_position = "topcenter"
    # )
) 

dev.off()


###################################################################
###################################################################

#Quartile plus 208 TFs based clustering minus 20 original factors

###################################################################
###################################################################

output_dir <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/misc/TFs_Annotation_file.txt", sep="\t")

# 2kb up and downstream heatmap:
input_file <- "/Users/suryachhetri/Dropbox/for_chris/quartile_motif_files_memechip/quartilebased_plus208_minus20_ideastyle_promoter_state_dist_for_heatmap.txt"

read_file <- fread(input_file, sep="\t", header=TRUE) 
read_file$Total_Frac <- NULL #delete last column
read_file$REST_STATE <- NULL

read_df <- as.data.frame(read_file)
data <- read_df[,2:ncol(read_df)] # Subtract last 4 columns
mat_data <- as.matrix(data)

# Rowwise scaling :As interested in binding variance across the cis-regulatory region
# for any given factor
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)

output_file_name <- file.path(output_dir, "Quartilebased_plus208_minus20_IDEAS_state_heatmap_clustering.pdf")       
#output_file_name <- file.path(output_dir, "Promoter_state_heatmap_figure_minus_rest.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS States",
    row_title="Quartile based Factors plus 208 minus 20",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.2),
    column_names_gp = gpar(fontsize = 6)
    # km=4,
    ## For heatmap and annotation legends (add new-legends)
    # heatmap_legend_param = list(
    # title_gp = gpar(fontsize = 10), 
    # labels_gp = gpar(fontsize = 6),
    # legend_direction = "horizontal",
    # legend_width = unit(3, "cm"), title_position = "topcenter"
    # )
) 

dev.off()


