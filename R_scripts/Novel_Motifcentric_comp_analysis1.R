library(data.table)
library(dplyr)
library(ggplot2)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
tomfile <- "Novel_motifs_tomtom_fraction_table.txt"
cisbp_tomfile <- "Novel_motifs_cisbp_tomtom_fraction_table.txt"

tomtom_df = fread(file.path(output_dir, cisbp_tomfile), sep="\t")
names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_from_cisBP.pdf")       
pdf(output_file_name)

x_cutoff <- 0.05
y_cutoff <- last(tomtom_df[tomtom_df$tom_Evalue <= x_cutoff]$motif_fraction)
# [1] 0.84
ggplot(tomtom_df, aes(x=tom_Evalue, y=motif_fraction, color="Cum Motif Fraction")) +
        geom_point(size=1) + 
        geom_rug(color="darkgray") + theme_bw() +
        geom_vline(xintercept=x_cutoff,linetype=2) +
        geom_hline(yintercept=y_cutoff,linetype=2) +
        labs(x="Tomtom similarity Evalue", y="Fraction of Motif", color='') + # xlabel, ylabel and legend title
        scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        scale_y_continuous(breaks=seq(0.00,1,0.2)) +
        ggtitle("CIS-BP Database") +
        theme(plot.title = element_text(hjust = 0.5, face="bold"))
        # ylim(0,1)
dev.off()

########################################
# For jaspar2018:

jaspar_tomfile <- "Novel_motifs_jaspar2018_tomtom_fraction_table.txt"
tomtom_df = fread(file.path(output_dir, jaspar_tomfile), sep="\t")
names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_from_Jaspar2018.pdf")       
pdf(output_file_name)

x_cutoff <- 0.05
y_cutoff <- last(tomtom_df[tomtom_df$tom_Evalue <= x_cutoff]$motif_fraction)
# [1] 0.82
ggplot(tomtom_df, aes(x=tom_Evalue, y=motif_fraction, color="Cum Motif Fraction")) +
        geom_point(size=1) + 
        geom_rug(color="darkgray") + theme_bw() +
        geom_vline(xintercept=x_cutoff,linetype=2) +
        geom_hline(yintercept=y_cutoff,linetype=2) +
        labs(x="Tomtom similarity Evalue", y="Fraction of Motif", color='') + # xlabel, ylabel and legend title
        scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        scale_y_continuous(breaks=seq(0.00,1,0.2)) +
        ggtitle("JASPAR 2018 Database") +
        theme(plot.title = element_text(hjust = 0.5, face="bold"))
        # ylim(0,1)
dev.off()
 
#######################################
# Jaspar 2016:

jaspar_tomfile <- "Novel_motifs_tomtom_fraction_table.txt"
tomtom_df = fread(file.path(output_dir, jaspar_tomfile), sep="\t")
names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_from_Jaspar2016.pdf")       
pdf(output_file_name)

x_cutoff <- 0.05
y_cutoff <- last(tomtom_df[tomtom_df$tom_Evalue <= x_cutoff]$motif_fraction)
# [1] 0.81
ggplot(tomtom_df, aes(x=tom_Evalue, y=motif_fraction, color="Cum Motif Fraction")) +
        geom_point(size=1) + 
        geom_rug(color="darkgray") + theme_bw() +
        geom_vline(xintercept=x_cutoff,linetype=2) +
        geom_hline(yintercept=y_cutoff,linetype=2) +
        labs(x="Tomtom similarity Evalue", y="Fraction of Motif", color='') + # xlabel, ylabel and legend title
        scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        scale_y_continuous(breaks=seq(0.00,1,0.2)) +
        ggtitle("JASPAR 2016 Database") +
        theme(plot.title = element_text(hjust = 0.5, face="bold"))
        # ylim(0,1)
dev.off()

##########################################
# CisBP motif scan with identifier (Scatter plot):

cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df.txt"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top5.txt"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top10.txt"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top20.txt"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top50.txt"

tomtom_df = fread(file.path(output_dir, cisbp_scan_tomfile), sep="\t")
names(tomtom_df)[names(tomtom_df) == 'tom_E-value_y'] <- "tom_E_value_final"
names(tomtom_df)[names(tomtom_df) == 'cmo_E-value'] <- "cmo_E_value"
#names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs.pdf")       
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_top5tomtomHitsCompare.pdf")       
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_top10tomtomHitsCompare.pdf")       
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_top20tomtomHitsCompare.pdf")       
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_top50tomtomHitsCompare.pdf")       

pdf(output_file_name)
#tomtom_df$anno <- factor(tomtom_df$anno, levels=c("mismatch","match","no_match"))
ggplot(tomtom_df, aes(y=as.factor(tf_motif_id), x=tom_E_value_final, label = tf_motif_id )) +
        geom_point(aes(color=anno), size=0.9) + geom_text(check_overlap = TRUE, size=1.5) +
        scale_color_manual(values=c("green","orange","red")) +
        #geom_rug(color="darkgray", sides="b") +
        geom_vline(xintercept=x_cutoff,linetype=2) +
        #geom_hline(yintercept=y_cutoff,linetype=2) +
        labs(x="Tomtom similarity Evalue", y="DNA associated factors", color='') + # xlabel, ylabel and legend title
        scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        #scale_y_continuous(breaks=seq(0.00,1,0.2)) +
        ggtitle("CIS-BP Database") + theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, face="bold"),
	    	  axis.text.y = element_text(size=0),
			 legend.title=element_text(face="bold")
		)
		
#		guides(fill=guide_legend(title="Annotation")) + theme_bw()
	        # ylim(0,1)
dev.off()

# Centrimo value (off peak distance):
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_offpeak_dist.pdf")       
pdf(output_file_name)
# Excluding "no_match"
tomtom_df_subset = tomtom_df[!tomtom_df$anno == "no_match",]
match_subset = tomtom_df[tomtom_df$anno == "match",]
mismatch_subset = tomtom_df[tomtom_df$anno == "mismatch",]
ks_test = ks.test(match_subset$cmo_bin_location,mismatch_subset$cmo_bin_location)
ggplot(tomtom_df_subset, aes(y=as.factor(tf_motif_id), x=cmo_bin_location, label = tf_motif_id )) +
        geom_point(aes(color=anno), size=0.9) + geom_text(check_overlap = TRUE, size=1.5) +
        scale_color_manual(values=c("green","orange","red")) +
        #geom_rug(color="darkgray", sides="b") +
        geom_vline(xintercept=0,linetype=2) +
        #geom_hline(yintercept=y_cutoff,linetype=2) +
        labs(x="Centrimo Offpeak Distance", y="DNA associated factors", color='') + # xlabel, ylabel and legend title
        #scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        #scale_y_continuous(breaks=seq(0,50,5)) +
        ggtitle(paste("Kolmogorov-Smirnov Test (pval) :", round(ks_test$p.value,4))) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold"),
	    	  axis.text.y = element_text(size=0),
			 legend.title=element_text(face="bold")
		) + facet_wrap(~anno)
		
#		guides(fill=guide_legend(title="Annotation")) + theme_bw()
	        # ylim(0,1)
dev.off()

# Custom inhouse script Centrimo value (off peak distance):
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset_nearest.txt"
tomtom_df = fread(file.path(output_dir, cisbp_scan_tomfile), sep="\t")
names(tomtom_df)[names(tomtom_df) == 'tom_E-value_y'] <- "tom_E_value_final"
names(tomtom_df)[names(tomtom_df) == 'cmo_E-value'] <- "cmo_E_value"
#names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_offpeak_dist_nearest.pdf")       
pdf(output_file_name)
# Excluding "no_match"
tomtom_df_subset = tomtom_df[!tomtom_df$anno == "no_match",]
match_subset = tomtom_df[tomtom_df$anno == "match",]
mismatch_subset = tomtom_df[tomtom_df$anno == "mismatch",]
ks_test = ks.test(match_subset$adj_offset_dist,mismatch_subset$adj_offset_dist)
ggplot(tomtom_df_subset, aes(y=as.factor(tf_motif_id), x=adj_offset_dist, label = tf_motif_id )) +
        geom_point(aes(color=anno_new), size=0.9) + geom_text(check_overlap = TRUE, size=1.5) +
        scale_color_manual(values=c("green","orange","red")) +
        #geom_rug(color="darkgray", sides="b") +
        geom_vline(xintercept=0,linetype=2) +
        geom_vline(xintercept=14,linetype=2) +
        labs(x="Offpeak median Distance", y="DNA associated factors", color='') + # xlabel, ylabel and legend title
        #scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        #scale_y_continuous(breaks=seq(0,50,5)) +
        ggtitle(paste("Kolmogorov-Smirnov Test (pval) :", 3.84e-11)) + #ks_test$p.value
        theme_bw()+ theme(plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.y = element_text(size=0),
             legend.title=element_text(face="bold")
        ) + facet_wrap(~anno_new) 
        
#       guides(fill=guide_legend(title="Annotation")) + theme_bw()
            # ylim(0,1)
dev.off()


####################################################

# Motif and Factors relation with Prom/Enh/Prom-Enh

####################################################


### library("readxl") Read xls with IDEAS annotation table:
# ideas_table <- read_excel(file.path(output_dir, "IDEAS_reannotation_table.xlsx" ))
# ideas_table_df <- ideas_table %>% as.data.frame

motif_factor_df <- fread(file.path(output_dir, "Motif_factor_and_regulatory_region_association.txt"))
output_file_name <- file.path(output_dir, "Motif_factor_and_regulatory_region_association.pdf")       
pdf(output_file_name)

ggplot(motif_factor_df, aes(x=anno_new, fill=final_ano)) +
  geom_bar(width=0.5) + scale_y_continuous(breaks=seq(0,200,20)) +
  labs(x="", y="Count", fill="Annotation") #+
  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")

dev.off()

        

##########################################
# Motif similarity Heatmap data:

######################################
#### Combined heatmap:

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
output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"

# Load TF annoation file:
# anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
# anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))


# Motif1
read_file <- fread(file.path(output_dir, "Novel_motifs_combined_tomtom_cisbp_allscan_df_heatmap.txt"), sep="\t", header=TRUE)
# read_file_merged <- merge(anno_df, read_file, by.x="Target", by.y="tf_name", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file)
#read_df <- read_df[apply(read_df[,1:ncol(read_df)], MARGIN = 1, FUN = function(x) sd(x) != 0),]

# read_df <- data.frame(sapply(read_df[,1:ncol(read_df)], function(x) as.numeric(as.character(x))))
#read_df <- read_df %>% filter(Category=="DBF")
#data <- read_df[,-c(1,2:6,7,8,9)] # Subtract last 4 columns
#data <- data %>% select("Heterochrom_repressed", "Weak Enh", "Strong Enh", "Prom assoc") # rearrange cols
mat_data <- data.matrix(read_df[,1:ncol(read_df)])
mat_data[is.na(mat_data)] <- 2.0

# Which motifs has the highest fraction say so normalized across the rows:
# mat_data1 <- t(apply(mat_data1, 1, scale))


colnames(mat_data) <- colnames(read_df)
rownames(mat_data) <- read_df$tf_motif_id # read_df$tf_name_edit
rownames(mat_data)
colnames(mat_data)

output_file_name <- file.path(output_dir, "Combined_tomtom_cisbp_heatmap_original.pdf")       
pdf(output_file_name)

# ht1 = Heatmap(mat_data, name="CisBinding Fraction", 
#     column_title="Motif1",
#     row_title="Transcription Factors",
#     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     row_names_gp = gpar(fontsize = 1.4),
#     column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

# ht2 = Heatmap(mat_data2, name="CisBinding Fraction", 
#     column_title="Motif2",
#     row_title="Transcription Factors",
#     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     row_names_gp = gpar(fontsize = 1.4),
#     column_names_gp = gpar(fontsize = 6), km=4, cluster_columns = FALSE) 

# ht_anno = Heatmap(read_df$annotation, name="annotation", width=unit(5, "mm"), 
#     heatmap_legend_param = list(
#     title_gp = gpar(fontsize = 10, fontface = "bold"), 
#     labels_gp = gpar(fontsize = 6, fontface = "bold")
#     ), col = c("Ideas Prom" = "red", "Ideas Enh" = "orange", "Ideas CTCF" = "blueviolet", "Ideas multiple" = "burlywood", "Ideas Other" = "green"))


# ht1 + ht2 + ht_anno


ht1 = Heatmap(mat_data, name="Similarity Eval", 
		col = colorRamp2(c(0.025,0.05,1,1.1), c("green", "blue", "red","grey")),
		column_title="CIS-BP Motifs",
		row_title="ChIPSeq Motifs",
		row_title_gp = gpar(fontsize = 15, fontface = "bold"),
		column_title_gp = gpar(fontsize = 15, fontface = "bold"),
		row_names_gp = gpar(fontsize = 1),
		column_names_gp = gpar(fontsize = 1), 
		na_col = "orange") 

ht1


dev.off()


###############################

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(cluster)

# args <-  commandArgs(TRUE)
# input_file <- args[1]
# tf_name <- args[2]
# out_dir_name <- args[3]

# Output dir:
output_dir <- "~/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
read_file <- fread(file.path(output_dir, "Novel_motifs_combined_tomtom_cisbp_allscan_df_heatmap.txt"), sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)
#read_df[is.na(read_df)] <- -0.05

data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data) 


output_file_name <- paste0("~/Dropbox", "/", "hepg2_liver_57_sp_gene_heatmap.pdf")        
pdf(output_file_name)

summary(read_df) #max_val = 325
#pair_break =  c(c(-1,0.099), seq(0,1,length=9))
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.50)
heatmap.2(mat_data,
  main = "TF enrichment over highly liver specific & hepg2 expressed genes (-log2(pvalue))", 
  xlab = "Liver-HepG2 specific genes",
  ylab = "Transcription Factors",
  col=greenred(10),  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  #symbreaks = min(mat_data, na.rm=TRUE),
  na.color="grey"
#  breaks=c(c(-0.1,-0.02),seq(0,1,(1-0)/9))
  )    

dev.off()
