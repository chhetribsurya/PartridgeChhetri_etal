library(igraph)
library(gtools)
library(ggplot2)
library(ggnet) # library(GGally)		
library(qgraph)
library(dplyr)
library(data.table)

################################################

# Network Plots

################################################


library(igraph)
library(gtools)
library(ggplot2)
library(ggnet) # library(GGally)		
library(qgraph)
library(dplyr)
library(data.table)

# For all tfs
# net_input <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/old_dir/final_cobinding_peak_count_fullsorted.bed", sep='\t', header=T)
#net_input <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_co_binding_total/final_cobinding_peak_count_perc_overlap_sorted_all_tf.bed", sep='\t', header=T)
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



# For DNA binding factors
#net_input <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_co_binding_total/final_cobinding_peak_count_perc_overlap_sorted_dbf_tf.bed", sep='\t', header=T)
net_input <- fread("~/Dropbox/encode_3/tf_cobind_total/files/tf_cobind_data/final_cobinding_peak_count_perc_overlap_sorted_dbf_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
colnames(read_df)
read_df <- read_df[read_df$percent_overlap >=75,]

# """select only the TF with peak count more than 500"""
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
#E(net)[E(net)$percent_overlap >= 70 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 75 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 80 & E(net)$percent_overlap < 90]$color <- "green3"
E(net)[E(net)$percent_overlap >= 90]$color <- "red"

 
pdf("~/Dropbox/final_TF_cobinding_network_dbf_tf_test.pdf", width=6, height=7 ) #call the png writer
pdf("~/Dropbox/final_TF_cobinding_network_dbf_tf_500peaks_2.pdf", width=6, height=7 ) #call the png writer
#pdf("~/Dropbox/final_TF_cobinding_network_dbf_tf_1.pdf", width=6, height=7 ) #call the png writer
par(mai=c(0,0,1,0))     #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(net,				#the graph to be plotted
layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
main='TF Cobinding Network Analysis',	#specifies the title
vertex.label.dist=0.1,			#puts the name labels slightly off the dots
vertex.color="orange",			#puts the name labels slightly off the dots
vertex.frame.color='orange', 
#vertex.frame.color='blue', 		#the color of the border of the dots 
vertex.label.color='black',		#the color of the name labels
vertex.label.font=2,			#the font of the name labels
#vertex.label= ifelse(degree(net) > 8, V(net)$name, NA),
vertex.label=V(net)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
#vertex.label.cex=0.3,
vertex.label.cex= ifelse(degree(net) > 8, 0.35, 0.25),
edge.arrow.size=.3,		#specifies the size of the font of the labels. can also be made to vary
edge.curved=rep(0.5) 
)

dev.off() #dont forget to close the device



#For chromatin Regulators
#net_input <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_co_binding_total/final_cobinding_peak_count_perc_overlap_sorted_cr_tf.bed", sep='\t', header=T)
net_input <- fread("~/Dropbox/encode_3/tf_cobind_total/files/tf_cobind_data/final_cobinding_peak_count_perc_overlap_sorted_cr_tf.bed", sep='\t', header=T)
read_df <- as.data.frame(net_input)
colnames(read_df)
read_df <- read_df[read_df$percent_overlap >=60,]

# """select only the TF with peak count more than 500"""
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
E(net)[E(net)$percent_overlap >= 60 & E(net)$percent_overlap < 70]$color <- "purple3"
E(net)[E(net)$percent_overlap >= 70 & E(net)$percent_overlap < 80]$color <- "yellow3"
E(net)[E(net)$percent_overlap >= 80 & E(net)$percent_overlap < 90]$color <- "green3"
E(net)[E(net)$percent_overlap >= 90]$color <- "red"

 
pdf("~/Dropbox/final_TF_cobinding_network_cr_tf_test.pdf", width=6, height=7 ) #call the png writer
pdf("~/Dropbox/final_TF_cobinding_network_cr_tf_500peaks.pdf", width=6, height=7 ) #call the png writer
pdf("~/Dropbox/final_TF_cobinding_network_cr_tf_1.pdf", width=6, height=7 ) #call the png writer
par(mai=c(0,0,1,0))     #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
plot(net,				#the graph to be plotted
layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
main='TF Cobinding Network Analysis',	#specifies the title
vertex.label.dist=0.1,			#puts the name labels slightly off the dots
vertex.color="orange",			#puts the name labels slightly off the dots
vertex.frame.color='orange', 
#vertex.frame.color='blue', 		#the color of the border of the dots 
vertex.label.color='black',		#the color of the name labels
vertex.label.font=2,			#the font of the name labels
#vertex.label= ifelse(degree(net) > 8, V(net)$name, NA),
vertex.label=V(net)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
#vertex.label.cex=0.3,
vertex.label.cex= ifelse(degree(net) > 8, 0.5, 0.4),
edge.arrow.size=.3,		#specifies the size of the font of the labels. can also be made to vary
edge.curved=rep(0.5) 
)

dev.off() #dont forget to close the device






