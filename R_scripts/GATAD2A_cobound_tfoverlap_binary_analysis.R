
library(gtools)
library(gplots)
library(data.table)
library(dplyr)

# Output dir:
output_dir <- "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/motifs_compinfo_analysis/gatad2a_cobind_analysis"

# Read files:
read_file <- fread(file.path(output_dir, "merged_gatad2a_loci_for_coboundTF_overlap.txt"), sep="\t")

read_df <- as.data.frame(read_file)
data <- read_df[,3:(ncol(read_df))] # Subtract last 4 columns
mat_data <- as.matrix(data)

colnames(mat_data) <- colnames(data)
rownames(mat_data) <- read_df$Loci 

#rownames(mat_data)
colnames(mat_data)


# For one extra Sidecol bar:
output_file_name <- file.path(output_dir, "gatad2a_cobind_analysis_final.pdf")       
pdf(output_file_name)

grey_col <-  colorRampPalette("grey")
red_col <- colorRampPalette(c("indianred1"))
custom_col <- c(grey_col(1), red_col(1))

par(cex.main=0.7)
heatmap.2(mat_data,
  dendrogram = "col",
  Rowv = TRUE, 
  Colv = TRUE,
  #scale="none",
  main = "GATAD2A analysis", 
  xlab = "GATAD2A associated 6 DNA binding factors",
  ylab = "GATAD2A Loci",
  col=greenred(10), 
  #col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  labRow = FALSE, #cexRow = 0.0,
  cexCol = 0.8,
  na.color="grey"
  )    

dev.off()



#################################################

# library(data.table)
# library(dplyr)
# library(ComplexHeatmap)
# library(circlize)
# library(colorspace)
# library(GetoptLong)

# # Output dir:
# output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis/gatad2a_cobind_analysis"

# # Read files:
# read_file <- fread(file.path(output_dir, "merged_gatad2a_loci_for_coboundTF_overlap.txt"), sep="\t")
# #anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))

# read_df <- as.data.frame(read_file)
# data <- read_df[,2:(ncol(read_df))] # Subtract last 4 columns
# mat_data <- as.matrix(data)

# colnames(mat_data) <- colnames(data)
# rownames(mat_data) <- read_df$Loci 

# rownames(mat_data)
# colnames(mat_data)


# # For one extra Sidecol bar:
# output_file_name <- file.path(output_dir, "fimo_motiffraction_heatmap_kmeans_plusDBF.pdf")       
# pdf(output_file_name)

# ht1 = Heatmap(mat_data, name="Binding Fraction", 
#     column_title="GATAD2A co-associated 6 DNA binding factors", na_col = "orange",
#     row_title="GATAD2A Loci",
#     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     row_names_gp = gpar(fontsize = 0),
#     column_names_gp = gpar(fontsize = 6), cluster_rows = TRUE,
#     cluster_columns = TRUE)

# # ht2 = Heatmap(read_df$annotation, name="category", width=unit(5, "mm"), 
# #     heatmap_legend_param = list(
# #     title_gp = gpar(fontsize = 10, fontface = "bold"), 
# #     labels_gp = gpar(fontsize = 6, fontface = "bold")
# #     ))
# # ht1 + ht2 
# ht1
# dev.off()



#############################################
