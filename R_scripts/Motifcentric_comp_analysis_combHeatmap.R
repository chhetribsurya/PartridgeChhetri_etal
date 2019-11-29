library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)
library(reshape)


###########################
# Complex heatmap Library
###########################


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

# Output dir:
output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))

# Read files:
read_file <- fread(file.path(output_dir, "final_fimomotif_peakfraction_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="TF_NAME", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data <- as.matrix(data)
mat_data[is.na(mat_data)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$TF_NAME # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)


# For one extra Sidecol bar:
output_file_name <- file.path(output_dir, "fimo_motiffraction_heatmap_kmeans_plusDBF.pdf")       
pdf(output_file_name)

ht1 = Heatmap(mat_data, name="Binding Fraction", 
    column_title="Motifs 1-5", na_col = "orange",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 2),
    column_names_gp = gpar(fontsize = 6), cluster_rows = TRUE,
    cluster_columns = TRUE)

# ht2 = Heatmap(read_df$annotation, name="category", width=unit(5, "mm"), 
#     heatmap_legend_param = list(
#     title_gp = gpar(fontsize = 10, fontface = "bold"), 
#     labels_gp = gpar(fontsize = 6, fontface = "bold")
#     ))
# ht1 + ht2 
ht1
dev.off()


#############
# Plot 2:
#############

# Read files:
read_file <- fread(file.path(output_dir, "MOTIF1_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data1 <- as.matrix(data)
mat_data1[is.na(mat_data1)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data1 <- t(apply(mat_data1, 1, scale))
colnames(mat_data1) <- colnames(data)

rownames(mat_data1) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data1)
colnames(mat_data1)

# For multi annotated column (2 cols extra):
output_file_name <- file.path(output_dir, "MOTIF1_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)

ht1 = Heatmap(mat_data1, name="CisBinding Z-score", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4) 
ht1
dev.off()


read_file <- fread(file.path(output_dir, "MOTIF2_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data2 <- as.matrix(data)
mat_data2[is.na(mat_data2)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data2 <- t(apply(mat_data2, 1, scale))
colnames(mat_data2) <- colnames(data)

rownames(mat_data2) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data2)
colnames(mat_data2)


output_file_name <- file.path(output_dir, "MOTIF2_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)
ht2 = Heatmap(mat_data2, name="CisBinding Z-score", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4)  

ht2
dev.off()

read_file <- fread(file.path(output_dir, "MOTIF3_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data3 <- as.matrix(data)
mat_data3[is.na(mat_data3)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data3 <- t(apply(mat_data3, 1, scale))
colnames(mat_data3) <- colnames(data)

rownames(mat_data3) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data3)
colnames(mat_data3)

output_file_name <- file.path(output_dir, "MOTIF3_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)

ht3 = Heatmap(mat_data3, name="CisBinding Fraction", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4) 

ht3

dev.off()

read_file <- fread(file.path(output_dir, "MOTIF4_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data4 <- as.matrix(data)
mat_data4[is.na(mat_data4)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data4 <- t(apply(mat_data4, 1, scale))
colnames(mat_data4) <- colnames(data)

rownames(mat_data4) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data4)
colnames(mat_data4)

output_file_name <- file.path(output_dir, "MOTIF4_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)

ht3 = Heatmap(mat_data4, name="CisBinding Fraction", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4) 

ht3

dev.off()

read_file <- fread(file.path(output_dir, "MOTIF5_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns
mat_data5 <- as.matrix(data)
mat_data5[is.na(mat_data5)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data5 <- t(apply(mat_data5, 1, scale))
colnames(mat_data5) <- colnames(data)

rownames(mat_data5) <- read_df$tf_name # read_df$tf_name_edit
rownames(mat_data5)
colnames(mat_data5)

output_file_name <- file.path(output_dir, "MOTIF5_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)

ht3 = Heatmap(mat_data5, name="CisBinding Fraction", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4) 

ht3

dev.off()


######################################
#### Combined heatmap:

# Motif1
read_file <- fread(file.path(output_dir, "MOTIF1_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data1 <- as.matrix(data)
mat_data1[is.na(mat_data1)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data1 <- t(apply(mat_data1, 1, scale))
colnames(mat_data1) <- colnames(data)

rownames(mat_data1) <- read_df$Target # read_df$tf_name_edit
rownames(mat_data1)
colnames(mat_data1)


# Motif2
read_file <- fread(file.path(output_dir, "MOTIF2_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data2 <- as.matrix(data)
mat_data2[is.na(mat_data2)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data2 <- t(apply(mat_data2, 1, scale))
colnames(mat_data2) <- colnames(data)

rownames(mat_data2) <- read_df$Target # read_df$tf_name_edit
rownames(mat_data2)
colnames(mat_data2)

# Motif3
read_file <- fread(file.path(output_dir, "MOTIF3_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data3 <- as.matrix(data)
mat_data3[is.na(mat_data3)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data3 <- t(apply(mat_data3, 1, scale))
colnames(mat_data3) <- colnames(data)

rownames(mat_data3) <- read_df$Target # read_df$tf_name_edit
rownames(mat_data3)
colnames(mat_data3)

# Motif4
read_file <- fread(file.path(output_dir, "MOTIF4_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data4 <- as.matrix(data)
mat_data4[is.na(mat_data4)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data4 <- t(apply(mat_data4, 1, scale))
colnames(mat_data4) <- colnames(data)

rownames(mat_data4) <- read_df$Target # read_df$tf_name_edit
rownames(mat_data4)
colnames(mat_data4)

# Motif5
read_file <- fread(file.path(output_dir, "MOTIF5_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.txt"), sep="\t", header=TRUE)
read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
read_df <- read_df %>% filter(Category=="DBF")
data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data5 <- as.matrix(data)
mat_data5[is.na(mat_data5)] <- 0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data5 <- t(apply(mat_data5, 1, scale))
colnames(mat_data5) <- colnames(data)

rownames(mat_data5) <- read_df$Target # read_df$tf_name_edit
rownames(mat_data5)
colnames(mat_data5)


output_file_name <- file.path(output_dir, "MOTIF_1234_final_tfcombined_fimocoords_with_ideascisreg_fractionbind_heatmapdata.pdf")       
pdf(output_file_name)

ht1 = Heatmap(mat_data1, name="CisBinding Fraction", 
    column_title="Motif1",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

ht2 = Heatmap(mat_data2, name="CisBinding Fraction", 
    column_title="Motif2",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 


ht3 = Heatmap(mat_data3, name="CisBinding Fraction", 
    column_title="Motif3",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 


ht4 = Heatmap(mat_data4, name="CisBinding Fraction", 
    column_title="Motif4",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 


ht5 = Heatmap(mat_data5, name="CisBinding Fraction", 
    column_title="Motif5",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

ht_anno = Heatmap(read_df$annotation, name="annotation", width=unit(5, "mm"), 
    heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 6, fontface = "bold")
    ), col = c("Ideas Prom" = "red", "Ideas Enh" = "orange", "Ideas CTCF" = "blueviolet", "Ideas multiple" = "burlywood", "Ideas Other" = "green"))


ht1 + ht2+ ht3 + ht4 + ht5+ ht_anno

dev.off()



output_file_name <- file.path(output_dir, "MOTIF_123_with_ideascisreg_fractionbind_heatmap_kmeans.pdf")       
pdf(output_file_name)

ht1 = Heatmap(mat_data1, name="CisBinding Fraction", 
    column_title="Motif1",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

ht2 = Heatmap(mat_data2, name="CisBinding Fraction", 
    column_title="Motif2",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 


ht3 = Heatmap(mat_data3, name="CisBinding Fraction", 
    column_title="Motif3",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

ht1 + ht2+ ht3 

dev.off()
