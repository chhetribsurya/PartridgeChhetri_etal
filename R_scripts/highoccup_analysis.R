library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)
library(gridExtra)
library(cowplot)

input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files/merged_d_100/final_tf_hotspots_sites_distribution_sorted_100bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)

tf_hotspot_count_df <-  tf_hotspot_df %>% 
					group_by(tf_counts, state_anno) %>%
					summarise(site_count=n()) %>% 
					data.frame

#zero_tf_df <- data.frame(tf_counts = 0, state_anno = "unbound_sites", site_count = 250179 )
zero_tf_df_prom <- data.frame(tf_counts = 0, state_anno = "prom_assoc", site_count = 23258 )
zero_tf_df_sEnh <- data.frame(tf_counts = 0, state_anno = "sEnh_assoc", site_count = 13037 )
zero_tf_df_wEnh <- data.frame(tf_counts = 0, state_anno = "wEnh_assoc", site_count = 213884 )

final_hotspots_df <- rbind(zero_tf_df_prom, zero_tf_df_sEnh, zero_tf_df_wEnh, tf_hotspot_count_df)
final_hotspots_df$log2_site_count <-  log2(final_hotspots_df$site_count)+1

final_hotspots_df %>% head(10)

write.table(final_hotspots_df, "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files/merged_d_100/final_tf_hotspots_sites_distribution_improved.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

prom_col <- colorRampPalette(c("red"))
strong_enh_col <- colorRampPalette(c("orange3"))
weak_enh_col <- colorRampPalette(c("yellow3"))
#zero_region <-  colorRampPalette("grey")

custom_col <- c(prom_col(1), strong_enh_col(1), weak_enh_col(1))
#custom_col <- c(zero_region(1), prom_col(1), strong_enh_col(1), weak_enh_col(1))

#read_file$hex_code <- custom_col

output_file_name <- paste0("~/Dropbox", "/", "Tf_hotspots_regulatory_region_distribution_merged_100bp_final_1.pdf")				
pdf(output_file_name)
# test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(final_hotspots_df, aes(x=tf_counts, y=log2_site_count, fill=state_anno)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across IDEAS HepG2 sites") + 
    theme(
    axis.text.y = element_text(size=8, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) +

	#theme(axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust = 1)) +
	#coord_flip() + 
	scale_fill_manual(values=custom_col )+
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() 

print(barplot)
dev.off()


### Stacked barplot for unmerged sites:

input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files/merged_d_1/final_tf_hotspots_sites_distribution_sorted_1bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)

tf_hotspot_count_df <-  tf_hotspot_df %>% 
					group_by(tf_counts, state_anno) %>%
					summarise(site_count=n()) %>% 
					data.frame


zero_tf_df <- data.frame(tf_counts = 0, state_anno = "unbound_sites", site_count = 343631 )
final_hotspots_df <- rbind(zero_tf_df, tf_hotspot_count_df) 
final_hotspots_df$log2_site_count <-  log2(final_hotspots_df$site_count)+1

final_hotspots_df %>% head(10)

write.table(final_hotspots_df, "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files/merged_d_1/final_tf_hotspots_sites_state_sorted.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

output_file_name <- paste0("~/Dropbox", "/", "final_tf_hotspots_regulatory_region_distribution_merged_1bp_1.pdf")				
pdf(output_file_name)
# test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(final_hotspots_df, aes(x=as.factor(tf_counts), y=log2_site_count, fill=state_anno)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across IDEAS HepG2 sites") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() 

print(barplot)
dev.off()


##############################
# Typical joint plots:


### For merged site, that is merging all the features within 100bp regions:
# input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_1/final_tf_hotspots_barplot_data.bed"
#input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


tf_hotspot_df <- tf_hotspot_df %>% select(tf_counts, total_sites, tf_counts_per_kb_segment)
tf_hotspot_df$log2_total_sites <-  log2(tf_hotspot_df$total_sites) + 1
tf_hotspot_df_1 <-  tf_hotspot_df %>% select(tf_counts, log2_total_sites) %>% mutate(Annotation="Site Counts")
tf_hotspot_df_2 <-  tf_hotspot_df %>% select(tf_counts, tf_counts_per_kb_segment) %>% mutate(Annotation="TF counts per kb")

### have same colnames to make rbind compatible; else rbind gets upset with different colnames:
names(tf_hotspot_df_1) <- c("tf_counts", "value", "Annotation")
names(tf_hotspot_df_2) <- c("tf_counts", "value", "Annotation")
combined_hotspots_df <- rbind(tf_hotspot_df_1, tf_hotspot_df_2)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

output_file_name <- paste0("~/Dropbox", "/", "final_tf_hotspots_and_tf_density_distribution_merged_100bp.pdf")				
pdf(output_file_name)

#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(combined_hotspots_df, aes(x=tf_counts, y=value, fill=Annotation)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts                            TF counts per kb segment") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across IDEAS HepG2 sites") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)

print(barplot)

dev.off()


### Corrected TF density plots:

### For merged site, that is merging all the features within 100bp regions:
# input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_1/final_tf_hotspots_barplot_data.bed"
#input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)
tf_hotspot_df$tf_counts_per_kb_segment_edited <- ((tf_hotspot_df$tf_counts * tf_hotspot_df$total_sites)/tf_hotspot_df$segment_size_sum)*1000

tf_hotspot_df$tf_counts_per_kb_segment_final <- with(tf_hotspot_df, tf_counts_per_kb_segment*total_sites)
tf_hotspot_df <- tf_hotspot_df %>% select(tf_counts, total_sites, tf_counts_per_kb_segment_final)
tf_hotspot_df$log2_total_sites <-  log2(tf_hotspot_df$total_sites) + 1
tf_hotspot_df_1 <-  tf_hotspot_df %>% select(tf_counts, log2_total_sites) %>% mutate(Annotation="Site Counts")
tf_hotspot_df_2 <-  tf_hotspot_df %>% select(tf_counts, tf_counts_per_kb_segment_final) %>% mutate(Annotation="TF counts per kb")

### have same colnames to make rbind compatible; else rbind gets upset with different colnames:
names(tf_hotspot_df_1) <- c("tf_counts_final", "value", "Annotation")
names(tf_hotspot_df_2) <- c("tf_counts_final", "value", "Annotation")
combined_hotspots_df <- rbind(tf_hotspot_df_1, tf_hotspot_df_2)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

#output_file_name <- paste0("~/Dropbox", "/", "final_tf_hotspots_and_tf_density_distribution_merged_100bp_corrected.pdf")				
output_file_name <- paste0("~/Dropbox/encode_3/geo_submission/paper/final_figures/supplementary", "/", "final_tf_hotspots_and_tf_density_plots.pdf")				
pdf(output_file_name)

#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(combined_hotspots_df, aes(x=tf_counts_final, y=value, fill=Annotation)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts                            TF counts per kb segment") + 
    theme_bw() + 
    ggtitle("TF bound sites distribution (across IDEAS HepG2 sites)") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() +  scale_x_reverse() + 
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)

print(barplot)

dev.off()

### Corrected TF site count and average size density plots:

### For merged site, that is merging all the features within 100bp regions:
# input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_1/final_tf_hotspots_barplot_data.bed"
#input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/final_tf_hotspots_barplot_data_100bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)
tf_hotspot_df$tf_sites_avg_size <- (tf_hotspot_df$segment_size_sum/tf_hotspot_df$total_sites)

#tf_hotspot_df$tf_counts_per_kb_segment_final <- with(tf_hotspot_df, tf_counts_per_kb_segment*total_sites)
tf_hotspot_df <- tf_hotspot_df %>% select(tf_counts, total_sites, tf_sites_avg_size)
tf_hotspot_df$log2_total_sites <-  log2(tf_hotspot_df$total_sites) + 1
tf_hotspot_df_1 <-  tf_hotspot_df %>% select(tf_counts, log2_total_sites) %>% mutate(Annotation="Site Counts")
tf_hotspot_df_2 <-  tf_hotspot_df %>% select(tf_counts, tf_sites_avg_size) %>% mutate(Annotation="Sites Avg. size(bp)")

### have same colnames to make rbind compatible; else rbind gets upset with different colnames:
names(tf_hotspot_df_1) <- c("tf_counts_final", "value", "Annotation")
names(tf_hotspot_df_2) <- c("tf_counts_final", "value", "Annotation")
combined_hotspots_df <- rbind(tf_hotspot_df_1, tf_hotspot_df_2)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

#output_file_name <- paste0("~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/plots_100bp_merged", "/", "final_tf_hotspots_and_avg_size_distribution_merged_100bp.pdf")				
output_file_name <- paste0("~/Dropbox/encode_3/geo_submission/paper/final_figures/supplementary", "/", "final_tf_hotspots_and_avg_size_distribution_merged_100bp.pdf")				
pdf(output_file_name)

#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

size_barplot <- ggplot(tf_hotspot_df_2, aes(x=tf_counts_final, y=value, fill=Annotation)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("Avg. size of the sites") + 
    theme_bw() + 
    #ggtitle("TF bound sites distribution (across IDEAS HepG2 sites)") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + scale_x_reverse() +
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)
#print(size_barplot)

sites_barplot <- ggplot(tf_hotspot_df_1, aes(x=tf_counts_final, y=value, fill=Annotation)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts") + 
    theme_bw() + 
    #ggtitle("TF bound sites distribution (across IDEAS HepG2 sites)") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + scale_x_reverse()+
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)

#print(sites_barplot)

plot_grid(sites_barplot, size_barplot, labels = c("A", "B"), align = "h")
#print(size_barplot)
#print(sites_barplot)
dev.off()


### For unmerged site, that is merged by 1bp overlap
#input_file <- "~/Dropbox/encode_3/tf_hotspots_total/files/merged_d_1/final_tf_hotspots_barplot_data.bed"
input_file <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files/merged_d_100/final_tf_hotspots_barplot_data_1bp.bed"
tf_hotspot_df <- fread(input_file, sep="\t", header= TRUE)


tf_hotspot_df <- tf_hotspot_df %>% select(tf_counts, total_sites, tf_counts_per_kb_segment)
tf_hotspot_df$log2_total_sites <-  log2(tf_hotspot_df$total_sites) + 1 
tf_hotspot_df_1 <-  tf_hotspot_df %>% select(tf_counts, log2_total_sites) %>% mutate(Annotation="Site Counts")
tf_hotspot_df_2 <-  tf_hotspot_df %>% select(tf_counts, tf_counts_per_kb_segment) %>% mutate(Annotation="TF counts per kb")

### have same colnames to make rbind compatible; else rbind gets upset with different colnames:
names(tf_hotspot_df_1) <- c("tf_counts", "value", "Annotation")
names(tf_hotspot_df_2) <- c("tf_counts", "value", "Annotation")
combined_hotspots_df <- rbind(tf_hotspot_df_1, tf_hotspot_df_2)

# change the order of factor in TF Category:
# meth_tf_df$tf_category <- factor(meth_tf_df$tf_category, levels = c("DBF", "CR/CF" ))

output_file_name <- paste0("~/Dropbox", "/", "final_tf_hotspots_and_tf_density_distribution_merged_1bp.pdf")				
pdf(output_file_name)

#test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))

barplot <- ggplot(combined_hotspots_df, aes(x=as.factor(tf_counts), y=value, fill=Annotation)) + 
	geom_bar(stat="identity") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    xlab("Transcription Factor Counts") + ylab("(Log2+1) site counts                            TF counts per kb segment") + 
    theme_bw() + 
    ggtitle("TF Hotspots distribution across IDEAS HepG2 sites") + 
    theme(
    axis.text.y = element_text(size=3, face="bold" ),
	plot.title=element_text(size=14, face="bold", hjust = 0.6),
	legend.title=element_text(face="bold")
	) + 
	coord_flip() + 
	guides(fill=guide_legend(title="Annotation")) + 
	scale_y_continuous() +
	facet_wrap(~Annotation)

print(barplot)

dev.off()


