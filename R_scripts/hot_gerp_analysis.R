library(dplyr)
library(ggplot2)


#input_dir <- "/Users/suryachhetri/Dropbox/encode_3/gerp_analysis/analysis_output"
input_dir <- "/Users/suryachhetri/Dropbox/encode_3/gerp_analysis/gerp_analysis_50bp_up_and_down/analysis_output"
input_file <- file.path(input_dir, "combined_chromwise_gerp_analysis_boxplot_data.txt")
output_dir <- file.path(input_dir, "plots")

## Provide the dir name(i.e sub dir) that you want to create under main dir:
#ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
if (!dir.exists(output_dir)){
	dir.create(output_dir)
} else {
	print("Dir already exists!")
}

gerp_df <- fread(input_file, sep="\t")

output_file_name <- file.path(output_dir, "tf_hotspot_sites_gerp_analysis_boxplot.pdf")			
pdf(output_file_name)

#desired_order <- c("1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
desired_order <- c("1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
gerp_df$binned_tf_count <- factor(gerp_df$binned_tf_count, levels=desired_order) #gerp_df$binned_tf_count %>% factor %>% levels

barplot <-ggplot(gerp_df, aes(x=binned_tf_count, y=mean_rs_score, fill=binned_tf_count)) + 
    geom_boxplot(width=0.6, fill = "grey80", outlier.colour = NA) + # outlier.shape = NA to remove the outliers 
	stat_boxplot( aes(binned_tf_count, mean_rs_score), 
    geom='errorbar', linetype=2, width=0.15, color="blue") +  
    stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
    xlab("Number of TFs bound") + ylab("GERP(RS Score)") + theme_bw()+ 
    ggtitle("GERP distribution on the TF hotspot sites") + 
    scale_y_continuous(limits=c(-2,2)) +
    theme(
    axis.text.x = element_text(size=10),
	plot.title=element_text(size=14,hjust = 0.6)
	)

print(barplot)
dev.off()

#######################################
#######################################

### Gerp element analysis per 100bp normalized:
input_dir <- "/Users/suryachhetri/Dropbox/encode_3/gerp_analysis/analysis_output"
input_file <- file.path(input_dir, "combined_gerp_element_analysis_boxplot_data_with_gerp_elem_overlap_fraction_final_edited.txt")
output_dir <- file.path(input_dir, "plots")

## Provide the dir name(i.e sub dir) that you want to create under main dir:
#ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
if (!dir.exists(output_dir)){
	dir.create(output_dir)
} else {
	print("Dir already exists!")
}

### Gerp element analysis per 100bp normalized with low bins:
gerp_df <- fread(input_file, sep="\t", header=TRUE)
gerp_df_final <- gerp_df %>% select(c("chrom", "start", "end", "anno", "TF_bound", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count"))
names(gerp_df_final) <- c("chrom", "start", "end", "anno_state", "tf_count", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count")

output_file_name <- file.path(output_dir, "Hotspot_sites_gerp_ELEMENT_boxplot_per100bp.pdf")			
pdf(output_file_name)

#desired_order <- c("1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
desired_order <- c("1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
gerp_df_final$binned_tf_count <- factor(gerp_df_final$binned_tf_count, levels=desired_order) #gerp_df_final$binned_tf_count %>% factor %>% levels

barplot <-ggplot(gerp_df_final, aes(x=binned_tf_count, y=as.numeric(fraction_overlap_per_100bp))) + 
    geom_boxplot(width=0.6, fill = "grey80", outlier.colour = NA) + # outlier.shape = NA to remove the outliers 
	stat_boxplot( aes(binned_tf_count, as.numeric(fraction_overlap_per_100bp)), 
    geom='errorbar', linetype=2, width=0.3, color="blue") +  
    #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Number of TFs bound") + ylab("Fraction overlap of GERP element per 100bp (Normalized)") + theme_bw()+ 
    ggtitle("GERP distribution on the TF hotspot sites") + 
    scale_y_continuous(limits=c(0,100)) +
    theme(
    axis.text.x = element_text(size=10),
	plot.title=element_text(size=14,hjust = 0.6)
	) 

#barplot <- barplot + geom_text(data = subset(gerp_df_final, meth_percent > 0.15), aes(x = tf_category, y =meth_percent, label = TF_name), 
#	hjust= -0.15, size = 1.5, check_overlap = TRUE)
print(barplot)

dev.off()


### Gerp element analysis per 100bp normalized with deep bins:
gerp_df <- fread(input_file, sep="\t", header=TRUE)
gerp_df_final <- gerp_df %>% select(c("chrom", "start", "end", "anno", "TF_bound", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count_zoomed"))
names(gerp_df_final) <- c("chrom", "start", "end", "anno_state", "tf_count", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count_zoomed")

output_file_name <- file.path(output_dir, "Hotspot_sites_gerp_ELEMENT_boxplot_per100bp_zoomed_bins.pdf")			
pdf(output_file_name)

#desired_order <- c("1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
desired_order <- c("1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100-149", "150+")
gerp_df_final$binned_tf_count_zoomed <- factor(gerp_df_final$binned_tf_count_zoomed, levels=desired_order) #gerp_df_final$binned_tf_count %>% factor %>% levels

barplot <-ggplot(gerp_df_final, aes(x=binned_tf_count_zoomed, y=as.numeric(fraction_overlap_per_100bp))) + 
    geom_boxplot(width=0.6, fill= "grey80", outlier.colour = NA) + # outlier.shape = NA to remove the outliers 
	stat_boxplot( aes(binned_tf_count_zoomed, as.numeric(fraction_overlap_per_100bp)), 
    geom='errorbar', linetype=2, width=0.3, color="blue") +  
    #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Number of TFs bound") + ylab("Fraction overlap of GERP element per 100bp (Normalized)") + theme_bw()+ 
    ggtitle("GERP distribution on the TF hotspot sites") + 
    scale_y_continuous(limits=c(0,100)) +
    theme(
    axis.text.x = element_text(size=10),
	plot.title=element_text(size=14,hjust = 0.6)
	) 

#barplot <- barplot + geom_text(data = subset(gerp_df_final, meth_percent > 0.15), aes(x = tf_category, y =meth_percent, label = TF_name), 
#	hjust= -0.15, size = 1.5, check_overlap = TRUE)

print(barplot)

dev.off()


#######################################
#######################################

### Gerp element analysis per 100bp normalized with deep bins and filled annotation for promoter and enhancers:
gerp_df <- fread(input_file, sep="\t", header=TRUE)
gerp_df_final <- gerp_df %>% select(c("chrom", "start", "end", "anno", "re_anno", "TF_bound", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count_zoomed"))
names(gerp_df_final) <- c("chrom", "start", "end", "anno_state", "re_anno_state", "tf_count", "fraction_overlap_per_100bp", "fraction_overlap_score_per_100bp", "binned_tf_count_zoomed")

output_file_name <- file.path(output_dir, "Hotspot_sites_gerp_ELEMENT_boxplot_per100bp_zoomed_bins_plot1.pdf")			
pdf(output_file_name)

#desired_order <- c("1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+")
desired_order <- c("1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100-149", "150+")
gerp_df_final$re_anno_state <- factor(gerp_df_final$re_anno_state, levels=c("Promoter_assoc", "Enhancer_assoc"))
gerp_df_final$binned_tf_count_zoomed <- factor(gerp_df_final$binned_tf_count_zoomed, levels=desired_order) #gerp_df_final$binned_tf_count %>% factor %>% levels

barplot <-ggplot(gerp_df_final, aes(x=binned_tf_count_zoomed, y=as.numeric(fraction_overlap_per_100bp), fill=re_anno_state)) + 
    geom_boxplot(width=0.6, outlier.colour = NA) + # outlier.shape = NA to remove the outliers 
	#stat_boxplot( aes(binned_tf_count_zoomed, as.numeric(fraction_overlap_per_100bp)), 
    #geom='errorbar', linetype=2, width=0.3, color="blue") +  
    #stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")+
    #geom_jitter(shape=16, position=position_jitter(0.004)) + 
    xlab("Number of TFs bound") + ylab("Fraction overlap of GERP element per 100bp (Normalized)") + theme_bw()+ 
    ggtitle("GERP distribution on the TF hotspot sites") + 
    scale_y_continuous(limits=c(0,100)) +
    theme(
    axis.text.x = element_text(size=10),
	plot.title=element_text(size=14,hjust = 0.6)
	) + 
	scale_fill_manual(values=c("indianred2", "orange")) +
	guides(fill=guide_legend(title="Annotation"))
	#facet_wrap (~re_anno_state, nrow=2)

print(barplot)

dev.off()
