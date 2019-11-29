library(igraph)
library(gtools)
library(ggplot2)
library(ggnet) # library(GGally)		
library(qgraph)
library(dplyr)
library(data.table)
library(stringr)

# df_10 = df[df["percent_overlap"] >=10]
# df_10 = df_10[df_10["total_peaks"] >=100]

""" For all tfs """
net_input <- fread("~/Dropbox/encode_3/tf_cobind_total/files/tf_cobind_data/final_cobinding_peak_count_perc_overlap_sorted_all_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
colnames(read_df)
read_df <- read_df[read_df$percent_overlap >=75,]

# select only the TF with peak count more than 500
read_df <- read_df[!read_df$total_peaks < 500,]
data <- read_df[,1:ncol(read_df)]

net<-graph.data.frame(data, directed=T)

V(net) # gives the list of vertices name which is single column only
E(net) # has all the attributes details with colnames; #attr: name (v/c), overlap_count (e/n), total_peaks (e/n)
el <- graph.edgelist(as.matrix(data[,2:3])) # gets the edge list of the columns you are interested in:
get.data.frame(net) # gets the list of all vertices and edges.
head(get.data.frame(net))

### set size of circle based on degree:
# V(net)$size # to see the list of size assigned to vertices or nodes
V(net)$size<-degree(net)/3 # here the size of the vertices is specified by the degree of the vertex, so that people supervising more have get proportionally bigger dots. Getting the right scale gets some playing around with the parameters of the scale function (from the 'base' package)

### color edges based on percent_overlap or weights:
#E(net)[E(net)$percent_overlap >= 60 & E(net)$percent_overlap < 70]$color <- "purple3"
E(net)[E(net)$percent_overlap >= 75 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 80 & E(net)$percent_overlap < 90]$color <- "green3"
E(net)[E(net)$percent_overlap >= 90]$color <- "red"

 
pdf("~/Dropbox/final_TF_cobinding_network_all_tf_peak_cutoffs_500_test.pdf", width=6, height=7 ) #call the png writer
pdf("~/Dropbox/final_TF_cobinding_network_all_tf_test.pdf", width=6, height=7 ) #call the png writer
par(mai=c(0,0,1,0))     #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(net,				#the graph to be plotted
layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
main='TF Cobinding Network Analysis',	#specifies the title
vertex.label.dist=0.1,
vertex.color="orange",			#puts the name labels slightly off the dots
vertex.frame.color='orange', 		#the color of the border of the dots 
vertex.label.color='black',		#the color of the name labels
vertex.label.font=2,			#the font of the name labels
#vertex.label= ifelse(degree(net) > 8, V(net)$name, NA),
vertex.label=V(net)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
#vertex.label.cex=0.3,
vertex.label.cex= ifelse(degree(net) > 10, 0.35, 0.25),
edge.arrow.size=.3,		#specifies the size of the font of the labels. can also be made to vary
edge.curved=rep(0.5) 
)

dev.off() #dont forget to close the device

###############################################

# Motif cobinding analysis:

""" For DNA binding factors """
net_input <- fread("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cobinding_analysis/dbf_tf/final_cobinding_motif_perc_overlap_sorted_dbf_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
read_df <- read_df %>% mutate(`tfmotif_id` = str_replace(`tfmotif_id`, "\\|MOTIF", ".M"))
read_df <- read_df %>% mutate(`cobind_tfmotif_id` = str_replace(`cobind_tfmotif_id`, "\\|MOTIF", ".M"))
colnames(read_df)

#read_df <- read_df[read_df$percent_overlap >=75,]
read_df <- read_df[read_df$percent_overlap >=70,]

# """select only the TF with peak count more than 500"""
# read_df <- read_df[!read_df$motif_sp_count < 500,]
data <- read_df[,2:ncol(read_df)]

net<-graph.data.frame(data, directed=T)

V(net) # gives the list of vertices name which is single column only
E(net) # has all the attributes details with colnames; #attr: name (v/c), overlap_count (e/n), motif_sp_count (e/n)
el <- graph.edgelist(as.matrix(data[,2:3])) # gets the edge list of the columns you are interested in:
get.data.frame(net) # gets the list of all vertices and edges.
head(get.data.frame(net))

### set size of circle based on degree:
# V(net)$size # to see the list of size assigned to vertices or nodes
V(net)$size<-degree(net)/3 # here the size of the vertices is specified by the degree of the vertex, so that people supervising more have get proportionally bigger dots. Getting the right scale gets some playing around with the parameters of the scale function (from the 'base' package)

### color edges based on percent_overlap or weights:
#E(net)[E(net)$percent_overlap >= 60 & E(net)$percent_overlap < 70]$color <- "purple3"
#E(net)[E(net)$percent_overlap >= 70 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 70 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 80 & E(net)$percent_overlap < 90]$color <- "green3"
E(net)[E(net)$percent_overlap >= 90]$color <- "red"

 
pdf("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cobinding_analysis/dbf_tf/final_motif_cobinding_network_dbf_test.pdf", width=6, height=7 ) #call the png writer

par(mai=c(0,0,1,0))     #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(net,				#the graph to be plotted
layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
main='TF Cobinding Network Analysis',	#specifies the title
vertex.label.dist=0.2,			#puts the name labels slightly off the dots
vertex.color="orange",			#puts the name labels slightly off the dots
vertex.frame.color='orange', 
#vertex.frame.color='blue', 		#the color of the border of the dots 
vertex.label.color='black',		#the color of the name labels
vertex.label.font=2,			#the font of the name labels
#vertex.label= ifelse(degree(net) > 8, V(net)$name, NA),
vertex.label=V(net)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
#vertex.label.cex=0.3,
vertex.label.cex= ifelse(degree(net) > 8, 0.25, 0.15),
edge.arrow.size=.2,		#specifies the size of the font of the labels. can also be made to vary
edge.curved=rep(0.5) 
)

dev.off() #dont forget to close the device

##########################################

# Motif co-occurence analysis

""" For DNA binding factors """
net_input <- fread("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cooccurrence_analysis/dbf_tf/final_cobinding_motif_perc_overlap_sorted_dbf_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
read_df <- read_df %>% mutate(`tfmotif_id` = str_replace(`tfmotif_id`, "\\|MOTIF", ".M"))
read_df <- read_df %>% mutate(`cobind_tfmotif_id` = str_replace(`cobind_tfmotif_id`, "\\|MOTIF", ".M"))
colnames(read_df)

#read_df <- read_df[read_df$percent_overlap >=75,]
read_df <- read_df[read_df$percent_overlap >=25,]

# """select only the TF with peak count more than 500"""
read_df <- read_df[!read_df$motif_sp_count < 100,]
data <- read_df[,2:ncol(read_df)]

net<-graph.data.frame(data, directed=T)

V(net) # gives the list of vertices name which is single column only
E(net) # has all the attributes details with colnames; #attr: name (v/c), overlap_count (e/n), motif_sp_count (e/n)
el <- graph.edgelist(as.matrix(data[,2:3])) # gets the edge list of the columns you are interested in:
get.data.frame(net) # gets the list of all vertices and edges.
head(get.data.frame(net))

### set size of circle based on degree:
# V(net)$size # to see the list of size assigned to vertices or nodes
V(net)$size<-degree(net)/3 # here the size of the vertices is specified by the degree of the vertex, so that people supervising more have get proportionally bigger dots. Getting the right scale gets some playing around with the parameters of the scale function (from the 'base' package)

### color edges based on percent_overlap or weights:
#E(net)[E(net)$percent_overlap >= 60 & E(net)$percent_overlap < 70]$color <- "purple3"
#E(net)[E(net)$percent_overlap >= 70 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 25 & E(net)$percent_overlap < 30]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 30 & E(net)$percent_overlap < 35]$color <- "green3"
E(net)[E(net)$percent_overlap >= 35]$color <- "red"

 
pdf("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cooccurrence_analysis/dbf_tf/final_motif_cooccurence_networkplot_dbf_100_25.pdf", width=6, height=7 ) #call the png writer

par(mai=c(0,0,1,0))     #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(net,				#the graph to be plotted
layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
main='TF Cobinding Network Analysis',	#specifies the title
vertex.label.dist=0.2,			#puts the name labels slightly off the dots
vertex.color="orange",			#puts the name labels slightly off the dots
vertex.frame.color='orange', 
#vertex.frame.color='blue', 		#the color of the border of the dots 
vertex.label.color='black',		#the color of the name labels
vertex.label.font=2,			#the font of the name labels
#vertex.label= ifelse(degree(net) > 8, V(net)$name, NA),
vertex.label=V(net)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
#vertex.label.cex=0.3,
vertex.label.cex= ifelse(degree(net) > 8, 0.30, 0.25),
edge.arrow.size=.2,		#specifies the size of the font of the labels. can also be made to vary
edge.curved=rep(0.5) 
)

dev.off() #dont forget to close the device


################################
# Misc:

net_input <- fread("~/Dropbox/for_genemodels/motifs_compinfo_analysis/cooccurrence_analysis/dbf_tf/final_cobinding_motif_perc_overlap_sorted_dbf_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
read_df <- read_df %>% mutate(`tfmotif_id` = str_replace(`tfmotif_id`, "\\|MOTIF", ".") %>% str_replace("_human", "") %>% str_replace("\\[FLAG\\]", ""))
read_df <- read_df %>% mutate(`cobind_tfmotif_id` = str_replace(`cobind_tfmotif_id`, "\\|MOTIF", ".") %>% str_replace("_human", "") %>% str_replace("\\[FLAG\\]", ""))
colnames(read_df)

# read_df <- read_df[read_df$percent_overlap >=75,]
# read_df <- read_df[read_df$percent_overlap >=25,]

# # """select only the TF with peak count more than 500"""
# read_df <- read_df[!read_df$motif_sp_count < 100,]
# data <- read_df[,2:ncol(read_df)]

# Filtered motif:
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
filt_df <- fread(file.path(output_dir, "Motif_factor_and_regulatory_region_association.txt"), sep="\t")

merged_df = merge(filt_df, read_df, by.x="tf_motif_id", by.y="tfmotif_id", all.x=TRUE)


########################################################################
# Motif cooccurrence heatmap:
# Multiple column based annotation heatmap - Complex Pheatmap params:

library(data.table)
library(dplyr)
library(pheatmap)
library(stringr)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
filt_motif_df <- fread(file.path(output_dir, "Motif_factor_and_regulatory_region_association.txt"), sep="\t")
heatmap_df <- fread(file.path(output_dir,"final_motif_cooccurrence_perc_overlap_heatmap_data.txt"), sep="\t", header=TRUE)
#heatmap_df <- fread(file.path(output_dir,"final_motif_cooccurrence_perc_overlap_heatmap_data_5percentmore.txt"), sep="\t", header=TRUE)
heatmap_df <- as.data.frame(heatmap_df)
mat <- heatmap_df[,2:ncol(heatmap_df)] %>% as.matrix

rownames(mat) <-  heatmap_df[,1]
colnames(mat)

# Filter motif 
final_filt_motif_df <- filt_motif_df[filt_motif_df$tf_motif_id %in% names(heatmap_df)]
final_filt_motif_df <- final_filt_motif_df %>% mutate(`final_ano` = str_replace(`final_ano`, "\\s", "_") %>% str_replace("-", "_") %>% str_replace("\\s", "_"))
#merge(colnames(mat))

motif_df <- final_filt_motif_df[,c("tf_motif_id", "anno_new", "final_ano", "Category")]

annotation_cols <- data.frame(annotation = motif_df$anno_new, cis_region = motif_df$final_ano)
annotation_rows <- data.frame(annotation = motif_df$anno_new, dbf_category = motif_df$Category)
#annotation_cols <- data.frame(annotation = motif_df$anno_new)

rownames(annotation_cols) <- motif_df$tf_motif_id
rownames(annotation_rows) <- motif_df$tf_motif_id

annotation_colors <- list(annotation=c(concordance="green", discordance="orange", no_match="red"), 
							cis_region=c(Promoter_like="red", Enhancer_like="orange3", Promoter_and_Enhancer_like="green", Heterochromatin_and_Repressed ="gray67", CTCF_Bound="purple3"),
							dbf_category=c(DBF="green", CR="black"))
#annotation_colors <- list(annotation=c(concordance="green", discordance="orange", no_match="red"))
col.pal <- RColorBrewer::brewer.pal(9, "Reds")

# library(grid)
# grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
# grid.gedit("col_annotation", gp = gpar(col="grey70"))

library(grid)
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
pheatmap(mat,
		scale="row",
		color = col.pal,
        annotation_col = annotation_cols, 
        annotation_row = annotation_rows,
        annotation_colors = annotation_colors,
        annotation_legend = TRUE,
        show_colnames = T, show_rownames = T, cluster_rows = T, 
        cluster_cols = T, legend = TRUE, 
        fontsize = 6,
        fontsize_row=1.2,
        fontsize_col = 1.2,
        cellheight = 1.2,
        cellwidth = 1,
        #filename = file.path(output_dir, "TEST1.pdf"),
        clustering_distance_rows = "euclidean", border_color = FALSE)

# gaps_row=c(10,20,30,40),
# gaps_col=c(3,6,9),
# cellheight = 6,
# cellwidth = 20,
# clustering_distance_rows = "correlation",#Pearson's
# clustering_method = "average",

###############################################################
###############################################################
# pheatmap for motif co-occurrence:
# Single column based annotation heatmap:

# Motif cooccurrence heatmap:
library(data.table)
library(dplyr)
library(pheatmap)
library(stringr)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
filt_motif_df <- fread(file.path(output_dir, "Motif_factor_and_regulatory_region_association.txt"), sep="\t")
heatmap_df <- fread(file.path(output_dir,"final_motif_cooccurrence_perc_overlap_heatmap_data.txt"), sep="\t", header=TRUE)
heatmap_df <- as.data.frame(heatmap_df)
mat <- heatmap_df[,2:ncol(heatmap_df)] %>% as.matrix

rownames(mat) <-  heatmap_df[,1]
colnames(mat)

# Filter motif 
final_filt_motif_df <- filt_motif_df[filt_motif_df$tf_motif_id %in% names(heatmap_df)]
final_filt_motif_df <- final_filt_motif_df %>% mutate(`final_ano` = str_replace(`final_ano`, "\\s", "_") %>% str_replace("-", "_") %>% str_replace("\\s", "_"))
motif_df <- final_filt_motif_df[,c("tf_motif_id", "anno_new", "final_ano", "Category")]

annotation_cols <- data.frame(Cis_Region = motif_df$final_ano)
annotation_rows <- data.frame(Cis_Region = motif_df$final_ano)

rownames(annotation_cols) <- motif_df$tf_motif_id
rownames(annotation_rows) <- motif_df$tf_motif_id

annotation_colors <- list(Cis_Region=c(Promoter_like="red", Enhancer_like="orange", Promoter_and_Enhancer_like="#619CFF", Heterochromatin_and_Repressed ="gray67", CTCF_Bound="purple3"))
col.pal <- RColorBrewer::brewer.pal(9, "Reds")
# col.pal <- RColorBrewer::brewer.pal(10, "RdBu")
###col.pal <- rev(col.pal) # for blue to red instead
###col.pal[6] <- "#619CFF" #col.pal <- c("blue", "red")
###col.pal <- col.pal[-7:-8]

# col.pal <- colorRampPalette(c("blue", "red"))(15)
# col.pal[8] <- "#619CFF"
# red_gradient <- rev(col.pal[9:length(col.pal)])
# blue_gradient <- rev(col.pal[1:7])
# mid <- "#619CFF"

# final_col.pal <- c(blue_gradient, mid, red_gradient)

library(grid)
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
pheatmap(mat,
		#scale="row",
		#color = final_col.pal,
		color = col.pal,
        annotation_col = annotation_cols, 
        #annotation_row = annotation_rows,
        annotation_colors = annotation_colors,
        annotation_legend = TRUE,
        show_colnames = T, show_rownames = T, cluster_rows = T, 
        cluster_cols = T, legend = TRUE, 
        fontsize = 6,
        fontsize_row=1.2,
        fontsize_col = 1.2,
        cellheight = 1.2,
        cellwidth = 1,
        #filename = file.path(output_dir, "Motif_cooccurrence_strength_heatmap_clustered.pdf"),
        clustering_distance_rows = "euclidean", border_color = FALSE,
        main = "Motif Co-occurence Heatmap")




