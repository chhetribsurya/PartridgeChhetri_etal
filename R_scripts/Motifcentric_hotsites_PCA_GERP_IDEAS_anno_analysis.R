
library(scatterplot3d) # load
library(data.table)
library(dplyr)
library(ggplot2)

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_30more_orig.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_30more_orig_wo_pca.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_50more_orig.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_50more_orig_wo_pca.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_10more_orig.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "tsne_analysis_hotmotif_sites_10more_orig_wo_pca.txt"), sep="\t", header=TRUE)

#######################
# 3-D tSNE plot
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_fig.pdf")
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_50more_figfinal.pdf")
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_50more_figfinal_wo_pca_new.pdf")
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_30more_figfinal_wo_pca.pdf")
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_10more_figfinal.pdf")
output_file_name = file.path("3D_tSNE_plot_Hotmotifsites_10more_figfinal_wo_pca.pdf")
pdf(file.path(output_dir, output_file_name))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
scatterplot3d(df[,1:3], main="3D t-SNE plot of Hotmotif sites", pch = 16, color=df$color_label, grid=FALSE, box=FALSE) 
addgrids3d(df[, 1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=c("Promoter","Strong Enhancer","Weak Enhancer","CTCF assoc","Other"),
col=c("red","orange", "yellow", "purple","green"), pch = 16, cex=1) 

dev.off()

############################
# 2-D tSNE plot

output_file_name <- "2D_tSNE_analysis_hotmotif_sites_fig.pdf"
output_file_name <- "2D_tSNE_analysis_hotmotif_sites_50more_figfinal_wo_pca_new.pdf"
output_file_name <- "2D_tSNE_analysis_hotmotif_sites_50more_figfinal.pdf"
output_file_name <- "2D_tSNE_analysis_hotmotif_sites_30more_figfinal_wo_pca.pdf"
output_file_name <- "2D_tSNE_analysis_hotmotif_sites_10more_figfinal.pdf"
output_file_name <- "2D_tSNE_analysis_hotmotif_sites_10more_figfinal_wo_pca.pdf"
pdf(file.path(output_dir, output_file_name), width=6.5, height=6)
df$ideas_state <-  factor(df$ideas_state, levels=c("Promoter","Strong Enhancer","Weak Enhancer","CTCF assoc","Other"))
plot <- ggplot(df, aes(x=tsne1, y=tsne2, color=ideas_state))+ 
		geom_point(size=1, alpha=0.7) +
		# geom_text(label = final_df["ideas_state"], size=6) +
		# geom_text(label=final_df["Label"], size=8)+
		ggtitle("t-SNE dimensions colored by Hotmotif site associated cis-region") +
		# theme_bw() +
		scale_x_continuous(name="x-tSNE_1") + scale_y_continuous(name = "y-tSNE_2" ) +
		# xlab("Principal Component 1") + ylab("Principal Component 2") +
		scale_color_manual(values=c("red", "orange", "yellow", "blueviolet", "green")) +
		guides(fill=guide_legend(title="IDEAS States")) +
		theme(
		axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		plot.title = element_text(size=12, face="bold"),
		legend.title = element_text(size=8, face="bold")
		# axis_text_y = element_text(size=1.3),
		# axis_text_y = element_text(size=1.3) 
		)
plot

dev.off()


########################
# PCA analysis of HotMotif sites:
########################

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
df = fread(file.path(output_dir, "PCA_analysis_hotmotif_sites_2.txt"), sep="\t", header=TRUE)
df = fread(file.path(output_dir, "PCA_analysis_hotmotif_sites_2_wo_scale.txt"), sep="\t", header=TRUE)

# 3-D PCA plots:
output_file_name = file.path("3D_PCA_plot_Hotmotifsites_fig.pdf")
output_file_name = file.path("3D_PCA_plot_Hotmotifsites_fig_wo_scale.pdf")
pdf(file.path(output_dir, output_file_name))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
scatterplot3d(df[,1:3], main="3D PCA plot of Hotmotif sites", pch = 16, color=df$color_label, grid=FALSE, box=FALSE) 
addgrids3d(df[, 1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=c("Promoter","Strong Enhancer","Weak Enhancer","CTCF assoc","Other"),
col=c("red","orange", "yellow", "purple","green"), pch = 16, cex=1) 

dev.off()

############################
# 2-D PCA plot
output_file_name <- "2D_PCA_analysis_hotmotif_sites_fig1.pdf"
output_file_name <- "2D_PCA_analysis_hotmotif_sites_fig_wo_scale.pdf"
pdf(file.path(output_dir, output_file_name), width=6.5, height=6)
# df$ideas_state <-  factor(df$ideas_state, levels=c("Other",  "CTCF assoc", "Promoter", "Weak Enhancer", "Strong Enhancer"))
df$ideas_state <-  factor(df$ideas_state, levels=c("Promoter","Strong Enhancer","Weak Enhancer","CTCF assoc","Other"))
plot <- ggplot(df, aes(x=Principal_component_1, y=Principal_component_2, color=ideas_state))+ 
		geom_point(size=1, alpha=0.7) +
		# geom_text(label = final_df["ideas_state"], size=6) +
		# geom_text(label=final_df["Label"], size=8)+
		ggtitle("PCA dimensions colored by Hotmotif site associated cis-region") +
		# theme_bw() +
		scale_x_continuous(name="Principal_component_1 (11.27%)") + scale_y_continuous(name = "Principal component 2 (3.85%)" ) +
		# xlab("Principal Component 1") + ylab("Principal Component 2") +
		# scale_color_manual(values=c("green", "blueviolet", "red", "yellow", "orange")) +
		scale_color_manual(values=c("red", "orange", "yellow", "blueviolet", "green")) +
		guides(fill=guide_legend(title="IDEAS States")) +
		theme(
		axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		plot.title = element_text(size=12, face="bold"),
		legend.title = element_text(size=8, face="bold")
		# axis_text_y = element_text(size=1.3),
		# axis_text_y = element_text(size=1.3) 
		)
plot

dev.off()


##################################
# Stacked barplot for Hotmotif sites annotation and their covered fraction:

output_file_name <- "Hotmotif_sites_annoation_barplot.pdf"
df = fread(file.path(output_dir, "Hotmotif_sites_annotation_and_fraction.txt"), sep="\t", header=TRUE)
pdf(file.path(output_dir, output_file_name), width=6.5, height=5.5)

custom_col <- c("red", "orange", "yellow", "purple", "green")
custom_col <- c("green", "purple", "yellow", "orange", "red")

#prom_col <- colorRampPalette(c("red"))
#strong_enh_col <- colorRampPalette(c("orange3"))
#weak_enh_col <- colorRampPalette(c("yellow3"))
#custom_col <- c(prom_col(1), strong_enh_col(1), weak_enh_col(1), "purple", "green")

df$ideas_state <-  factor(df$ideas_state, levels=c("Other", "CTCF assoc", "Weak Enhancer", "Strong Enhancer", "Promoter"))
barplot <- ggplot(df, aes(x=uniq_tfcount, y=element_perc, fill=ideas_state)) + 
	geom_bar(stat="identity", width = 1) +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Unique Factor Count") + ylab("Fraction of Cis-elements") + 
    theme_bw() + 
    ggtitle("Cis-element distribution at merged motif sites") + 
    theme(
    axis.text.y = element_text(size=8, face="bold" ),
	plot.title=element_text(size=12, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold"),
	# legend.position="top"
	) + #scale_x_reverse() + coord_flip() +
	scale_fill_manual(values=custom_col )+
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() 

print(barplot)
dev.off()


##################################
# Piechart : Top TF with recurring motif based on SVM score:

library(stringr)

output_file_name <- "Hotmotif_SVM_based_piechart.pdf"
output_file_name <- "Hotmotif_sites_problem3_piechart_data.pdf"
df = fread(file.path(output_dir, "Hotmotif_sites_problem3_piechart_data.txt"), sep="\t", header=TRUE)
pdf(file.path(output_dir, output_file_name), width=7, height=5.5)
df <- as.data.frame(df)

# df<-data.frame(table(TopHit))
col_names <- df[,1]
col_names <- col_names %>% str_replace("\\[FLAG\\]", "")
df[,6] <- col_names
# Cols<-rep("black",nrow(df))
# Cols[grep("FLAG",col_names)]<-"blue"
order_idx <- order(-df[,2])
num_counts <- df[order_idx,2]
total_labelcount <- 16
pie(num_counts, labels = c(col_names[1:total_labelcount], rep("",nrow(df)-total_labelcount)), font=1, main = "Top factors based on SVM score at HotMotif sites", cex=0.5, col = c("red","cornflowerblue"))

dev.off()







