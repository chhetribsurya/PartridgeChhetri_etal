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

output_dir <- "/Users/suryachhetri/Dropbox/paper_revision_3/lcp_hcp_prom_heatmap"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))

# 2kb up and downstream heatmap:
input_file <- "/Users/suryachhetri/Dropbox/paper_revision_3/lcp_hcp_prom_heatmap/ideastyle_promoter_state_dist_for_heatmap.txt"

# input_file <- "/home/surya/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap_zscore.txt"
# input_file <- "~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt"

read_file <- fread(input_file, sep="\t", header=TRUE)
read_file$Total_Frac <- NULL #delete last column
#read_file$REST_STATE <- NULL
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-5)] # Subtract last 4 columns
mat_data <- as.matrix(data)

# Rowwise scaling :As interested in binding variance across the cis-regulatory region
# for any given factor
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)

output_file_name <- file.path(output_dir, "Promoter_state_heatmap_figure.pdf")       
#output_file_name <- file.path(output_dir, "Promoter_state_heatmap_figure_minus_rest.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="Promoter States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
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


# For one extra Sidecol bar:
output_file_name <- file.path(output_dir, "Promoter_state_heatmap_plusDBF.pdf")       
#output_file_name <- file.path(output_dir, "Promoter_state_heatmap_plusDBF_minus_rest.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="Promoter States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6) #km=4
    ) + 
Heatmap(read_df$Category, name="category", width=unit(5, "mm"), 
    heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 6, fontface = "bold")
    ), col=c("DBF" = "green", "CR/CF" = "black")) 

# Heatmap(read_df$Category, name="category", width=unit(5, "mm"), 
#     heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
#     labels_gp = gpar(fontsize = 6)),
#     col = list(type2 = c("DBF" =  "green", "CR/CF" = "orange"))
#     )

dev.off()


##################################################
##################################################

# For multi annotated column (2 cols extra):
#output_file_name <- file.path(output_dir, "Promoter_state_piechart_heatmap_plusDBF_plusStates.pdf")       
#pdf(output_file_name)
#
#ht1 = Heatmap(mat_data, name="Binding Z-score", 
#    column_title="IDEAS Genomic States",
#    row_title="Transcription Factors",
#    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#    row_names_gp = gpar(fontsize = 1.4),
#    column_names_gp = gpar(fontsize = 6) #km=4
#    ) 
#
#ht2 = Heatmap(read_df$Category, name="category", width=unit(5, "mm"), 
#    heatmap_legend_param = list(
#    title_gp = gpar(fontsize = 10, fontface = "bold"), 
#    labels_gp = gpar(fontsize = 6, fontface = "bold")
#    ), col=c("DBF" = "green", "CR/CF" = "black")) 
#
#ht3 = Heatmap(read_df$annotation, name="annotation", width=unit(5, "mm"), 
#    heatmap_legend_param = list(
#    title_gp = gpar(fontsize = 10, fontface = "bold"), 
#    labels_gp = gpar(fontsize = 6, fontface = "bold")
#    ), col = c("Ideas Prom" = "red", "Ideas Enh" = "orange", "Ideas CTCF" = "blueviolet", "Ideas multiple" = "burlywood", "Ideas Other" = "green"))
#
#ht1 + ht3 + ht2
#
## ha_row = rowAnnotation(df = read_df, 
##     col = list(annotation = c("Prom assoc" = "red", "Enh assoc" = "orange", "CTCF assoc" = "blueviolet", "Multiple assoc" = "burlywood", "Others"="green")), width = unit(1, "cm"))
## ht_list = ht1 + ht2 + ha_row
## draw(ht_list)
#
#dev.off()
