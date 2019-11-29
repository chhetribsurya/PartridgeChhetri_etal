library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)
library(reshape)
# library(grid)

##################################
# HEATMAP.2 plots with clustering:
##################################

output_dir <- "~/Dropbox/for_genemodels"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

# 2kb up and downstream heatmap:
input_file <- "~/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap.txt"
# input_file <- "/home/surya/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap_zscore.txt"

read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] # Subtract last 4 columns

# Change everything to matrix dataset
mat_data <- as.matrix(data)
rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)


output_file_name <- file.path(output_dir, "TFs_peakfraction_5kb_up-downstream_tss.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("white", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  scale="column",
  main = "TFs binding distribution with genecode genes", 
  xlab = "Distance(bp) from the centre of TSS",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label, # matches with: read_df$tf_name to assign color
  key.xlab="Binding fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )   
 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()


### 3kb up and downstream heatmap:

### Load the TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")

### 1kb up and downstream heatmap:
input_file <- "~/Dropbox/encode_3/motifs_methplot_total/all_tf_files/all_tfs_motif_bin_methplot_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
#read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-4)] ## Subtract last 4 columns
mat_data <- as.matrix(data)

column_ranges <- as.numeric(colnames(mat_data))
column_select <- column_ranges[column_ranges >= -3000 & column_ranges <= 3000] 
mat_data <- mat_data[,colnames(mat_data) %in% column_select]

rownames(mat_data) <- read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)


output_file_name <- paste0(output_dir, "/", "final_motifs_binned_methplot_1kb_up&downstream.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(10)

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = FALSE,
  #scale="none",
  main = "Methylation distribution from the motif centres", 
  xlab = "Distance(bp) from the centre of motifs",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.35,
  na.color="grey",
  RowSideColors=read_df$label,
  key.xlab="Methylation fraction"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()


#########################################################################
# IDEAS piechart analysis-discrepancy impact-ratio and heatmap clustering
#########################################################################

#require("heatmap.plus")

output_dir <- "~/Dropbox/for_genemodels"

# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))
# 2kb up and downstream heatmap:
input_file <- "~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt"
# input_file <- "/home/surya/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap_zscore.txt"

read_file <- fread(input_file, sep="\t", header=TRUE)
read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged[,c("tf_name", "Category")]
read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file_merged)
data <- read_df[,2:(ncol(read_df)-5)] # Subtract last 4 columns

# Change everything to matrix dataset
mat_data <- as.matrix(data)
rownames(mat_data) <- read_df$tf_name # read_df$tf_name_edit

colnames(mat_data)
rownames(mat_data)


output_file_name <- file.path(output_dir, "Ideas_state_piechart_heatmap1.pdf")       
pdf(output_file_name)

#custom_col <- c(green_col(1), red_col1(1))
#custom_col <- c(grey_col(1), red_col(1))
custom_col <- colorRampPalette(c("blue", "yellow", "red"))(100)
#custom_col <- colorRampPalette(c("black", "white", "red"))(100)
myCols <- cbind(read_df$label2, read_df$label)
summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.80)
test <- heatmap.2(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = TRUE,
  scale="row",
  main = 'Different Chromatin "Flavors" for different TFs', 
  xlab = "IDEAS Genomic States",
  ylab = "Transcription Factors",
  #col=bluered(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #na.color="grey",
  RowSideColors=read_df$label, # matches with: read_df$tf_name to assign color
  key.xlab="Binding fraction",
  keysize=1 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          #key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          #lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )   
 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()




test <- heatmap.plus(mat_data,
  # dendrogram = "none",
  # Rowv = TRUE, 
  Colv = TRUE,
  scale="row",
  main = 'Different Chromatin "Flavors" for different TFs', 
  xlab = "IDEAS Genomic States",
  ylab = "Transcription Factors",
  #col=bluered(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #na.color="grey",
  RowSideColors=myCols, # matches with: read_df$tf_name to assign color
  key.xlab="Binding fraction",
  keysize=1 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          #key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          #lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1)
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )   
 
legend("topright", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
fill=c("red","orange","blueviolet","burlywood","green"), border=TRUE, bty="n", y.intersp = 0.9, cex=0.7) 

test
dev.off()


###########################
# Complex heatmap Library
###########################


library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)


# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))
# 2kb up and downstream heatmap:
input_file <- "~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt"
# input_file <- "/home/surya/Dropbox/for_genemodels/all_tf_tss_intersect_for_heatmap_zscore.txt"

read_file <- fread(input_file, sep="\t", header=TRUE)
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

output_file_name <- file.path(output_dir, "Ideas_state_piechart_heatmap_figure.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS Genomic States",
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
output_file_name <- file.path(output_dir, "Ideas_state_piechart_heatmap_plusDBF.pdf")       
pdf(output_file_name)

Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS Genomic States",
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
    )) 

# Heatmap(read_df$Category, name="category", width=unit(5, "mm"), 
#     heatmap_legend_param = list(title_gp = gpar(fontsize = 10), 
#     labels_gp = gpar(fontsize = 6)),
#     col = list(type2 = c("DBF" =  "green", "CR/CF" = "orange"))
#     )

dev.off()

# For multi annotated column (2 cols extra):
output_file_name <- file.path(output_dir, "Ideas_state_piechart_heatmap_plusDBF_plusStates.pdf")       
pdf(output_file_name)

ht1 = Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6) #km=4
    ) 

ht2 = Heatmap(read_df$Category, name="category", width=unit(5, "mm"), 
    heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 6, fontface = "bold")
    ), col=c("DBF" = "green", "CR/CF" = "black")) 

ht3 = Heatmap(read_df$annotation, name="annotation", width=unit(5, "mm"), 
    heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 6, fontface = "bold")
    ), col = c("Ideas Prom" = "red", "Ideas Enh" = "orange", "Ideas CTCF" = "blueviolet", "Ideas multiple" = "burlywood", "Ideas Other" = "green"))

ht1 + ht3 + ht2

dev.off()
