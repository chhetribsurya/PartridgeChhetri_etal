##########################################
# Fox Motif cobind Heatmap in R :
#########################################


library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

# Output dir:
output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis/fox_meme_motifs/occupancy_analysis"

# Read file:
read_file <- fread(file.path(output_dir, "Factors_occupancy_with_motif_details_37factors_fin.txt"), sep="\t", header=TRUE)

# Insert order for splitting the heatmap later:
read_file1 <- read_file[1:29,]
read_file2 <- read_file[30:37,]
read_file1$annotation = "Primary motif (+)" 
read_file2$annotation = "Primary motif (-)" 
final_read_file <- rbind(read_file1, read_file2)


read_df <- as.data.frame(read_file)
#read_df <- read_df %>% select("tf_name", "Self only", "Both motif", "No motif", "FOX only")

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
# mat_data1 <- t(apply(mat_data1, 1, scale))

rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)


output_file_name <- file.path(output_dir, "Factors_occupancy_with_motif_details_37factors_heatmap.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ht1 <- Heatmap(mat_data, name="Occupancy Fraction",
        #cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        col=color,
        column_title="Motif Status",
        row_title="Associated Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 5.5),
        column_names_gp = gpar(fontsize = 10), 
        na_col = "lightgray"
        #split = final_read_file$category) 
        )


ht_anno = Heatmap(final_read_file$annotation, name="", width=unit(5, "mm"), 
    heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 6, fontface = "bold")
    ), col = c("Primary motif (+)" = "lightgreen", "Primary motif (-)" = "grey"))

ht1 + ht_anno

dev.off()

##########################################
# Same Heatmap with box plot added though:
##########################################

read_file <- fread(file.path(output_dir, "FOXome_occupancy_with_follow_factors.txt"), sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
#read_df <- read_df %>% select("tf_name", "Self only", "Both motif", "No motif", "FOX only")

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
# mat_data1 <- t(apply(mat_data1, 1, scale))

rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)

output_file_name <- file.path(output_dir, "FOXome_occupancy_with_follow_factors_heatmap.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data, axis = TRUE, ylim=c(0,max(mat_data)), border=FALSE, gp = gpar(fill = "red")))
ht1 <- Heatmap(mat_data, name="Cooccupancy Fraction",
        #cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        col=color,
        column_title="Associated Factors",
        row_title="Forkhead Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 5), 
        na_col = "lightgray",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(3.8, "cm")) 

ht1

dev.off()


# For Barbara - Peak ordered heatmap:
read_file <- fread(file.path(output_dir, "FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
read_df <- read_df %>% filter(V1 != "tf_peakcount")
read_df <- read_df %>% rename("tf_name" = "V1")

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data1 <- t(apply(mat_data, 1, scale))
rownames(mat_data1) <- read_df$tf_name
colnames(mat_data1) <- colnames(read_df)[2:ncol(read_df)]
colnames(mat_data1)
rownames(mat_data1)

# For boxplot mat_data:
rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)

output_file_name <- file.path(output_dir, "FOXome_occupancy_with_follow_factors_heatmap_zscored.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data, axis = TRUE, ylim=c(0,max(mat_data)), border=FALSE, gp = gpar(fill = "red")))
ht1 <- Heatmap(mat_data1, name="Cooccupancy Fraction",
        #cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        # col=color,
        column_title="Associated Factors",
        row_title="Forkhead Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 5), 
        na_col = "lightgray",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(3.8, "cm")) 

ht1

dev.off()


# For mean_normalized - Peak ordered heatmap:
read_file <- fread(file.path(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
read_df <- read_df %>% filter(V1 != "tf_peakcount")
read_df <- read_df %>% rename("tf_name" = "V1")

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data1 <- t(apply(mat_data, 1, scale))
rownames(mat_data1) <- read_df$tf_name
colnames(mat_data1) <- colnames(read_df)[2:ncol(read_df)]
colnames(mat_data1)
rownames(mat_data1)

# For boxplot mat_data:
rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)

output_file_name <- file.path(output_dir, "mean_norm_FOXome_occupancy_with_follow_factors_heatmap_zscored.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data, axis = TRUE, ylim=c(0,max(mat_data)), border=FALSE, gp = gpar(fill = "red")))
ht1 <- Heatmap(mat_data1, name="Mean Norm Fraction",
        #cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        # col=color,
        column_title="Associated Factors",
        row_title="Forkhead Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 5), 
        na_col = "lightgray",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(3.8, "cm")) 

ht1

dev.off()


# For overlap size (ratio tendency to overlap) normalized - Peak ordered heatmap:
read_file <- fread(file.path(output_dir, "norm_FOXome_occupancy_with_follow_factors_peakOrdered.txt"), sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
read_df <- read_df %>% filter(V1 != "tf_peakcount")
read_df <- read_df %>% rename("tf_name" = "V1")

mat_data <- data.matrix(read_df[,2:ncol(read_df)])
# mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
mat_data1 <- t(apply(mat_data, 1, scale))
rownames(mat_data1) <- read_df$tf_name
colnames(mat_data1) <- colnames(read_df)[2:ncol(read_df)]
colnames(mat_data1)
rownames(mat_data1)

# For boxplot mat_data:
rownames(mat_data) <- read_df$tf_name
colnames(mat_data) 
rownames(mat_data)

output_file_name <- file.path(output_dir, "norm_FOXome_occupancy_with_follow_factors_unclust_peakOrdered.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Reds")

# Not Z-Scored though:
ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data, axis = TRUE, ylim=c(0,max(mat_data)), border=FALSE, gp = gpar(fill = "red")))
ht1 <- Heatmap(mat_data, name="Normalized Fraction",
        cluster_col=FALSE, 
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        col=color,
        column_title="Associated Factors",
        row_title="Forkhead Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 5), 
        na_col = "lightgray",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(3.8, "cm")) 

ht1

dev.off()
