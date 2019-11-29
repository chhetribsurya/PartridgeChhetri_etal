##########################################################
# 
# Motif cobinding fraction + Motif offset distance heatmap
#
##########################################################

# Read file:
read_file2 <- fread(file.path(output_dir, "Factors_occupancy_with_fox_motif_offsetdist.txt"), sep="\t", header=TRUE)

read_df2 <- as.data.frame(read_file2)

mat_data2 <- data.matrix(read_df2[,2:ncol(read_df2)])

rownames(mat_data2) <- read_df2$tf_name
colnames(mat_data2) 
rownames(mat_data2)

output_file_name <- file.path(output_dir, "Factors_cooccupancy_heatmap_with_fox_motifs_offsetdist_and_boxplot.pdf")       
pdf(output_file_name)

library(RColorBrewer)
color <- brewer.pal(9, "Blues")
color <- brewer.pal(7, "Spectral")
color <- colorRamp2(c(0,5,10,20,30,40,50,60,70,80), c("white", "darkgrey", "black", brewer.pal(7, "Spectral")))

ha_boxplot = HeatmapAnnotation(boxplot = anno_boxplot(mat_data2, axis = TRUE, ylim=c(0,max(mat_data2)), border=FALSE, gp = gpar(fill = "lightblue")))
ht2 <- Heatmap(mat_data2, name=" Offset Distance ", 
        cluster_col=FALSE,
        #col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
        #col = colorRamp2(c(0,1), c("blue", "red")),
        #col = colorRamp2(c(0,20,90), c("red", "lightblue", "steelblue")),
        #col = colorRamp2(c(0,20,90), brewer.pal(3, "Blues")),
        #col = colorRamp2(c(0,20,90), brewer.pal(3, "RdBu")), # colorRamp2(c(0,20,90)) - represents breaks
        col=color,
        column_title="FOX Motifs",
        row_title="Associated Factors",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_names_gp = gpar(fontsize = 4.5),
        column_names_gp = gpar(fontsize = 9), 
        na_col = "orange",
        bottom_annotation = ha_boxplot,
        bottom_annotation_height = unit(2.8, "cm")) 

ht2

dev.off()
