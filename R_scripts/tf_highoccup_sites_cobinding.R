library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

### Barplot for the TF covering the most no. of sites on TF hotspot genomic sites:
total_tf_hotspot_sites <- 12927
total_tf_hotspot_sites_enh <- 6227
total_tf_hotspot_sites_prom <- 6700

tf_hotspot_file <- fread("~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/tf_hotspots_cobinding_analysis_with_tf_counts.txt", sep="\t")
tf_hotspot_file_sorted <- tf_hotspot_file[order(tf_hotspot_file$count, decreasing=TRUE)]
tf_hotspot_file_sorted$sites_percent <- (tf_hotspot_file_sorted$count/total_tf_hotspot_sites)*100
tf_hotspot_file_sorted <- tf_hotspot_file_sorted %>% select(tf_name, count, sites_percent)
top_25_perc_whole_genome_assoc <- as.numeric(quantile(tf_hotspot_file_sorted$sites_percent, 0.75)) ## or, use [[]] to extract the value


tf_hotspot_region_based_file <- fread("~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/tf_hotspots_cobinding_analysis_with_region_based_tf_counts.txt", sep="\t")
tf_hotspot_region_based_file$sites_percent <- (tf_hotspot_region_based_file$count/total_tf_hotspot_sites)*100
tf_hotspot_region_based_file <- tf_hotspot_region_based_file %>% select(tf_name, state_anno, count, sites_percent)

### Filtering the enhancer associated region:
tf_hotspot_enh_associated <- tf_hotspot_region_based_file[tf_hotspot_region_based_file$state_anno == "sEnh_assoc"| tf_hotspot_region_based_file$state_anno == "wEnh_assoc"]
#tf_hotspot_enh_associated_df <- dcast(tf_hotspot_enh_associated, tf_name~state_anno, value.var="sites_percent")
tf_hotspot_enh_associated_df <- dcast(tf_hotspot_enh_associated, tf_name~state_anno, value.var="count")
tf_hotspot_enh_associated_df[is.na(tf_hotspot_enh_associated_df)] <- 0
tf_hotspot_enh_associated_df$total_enh_sites <- tf_hotspot_enh_associated_df$sEnh_assoc + tf_hotspot_enh_associated_df$wEnh_assoc
tf_hotspot_enh_associated_df$sites_percent <- (tf_hotspot_enh_associated_df$total_enh_sites/total_tf_hotspot_sites_enh)*100
tf_hotspot_enh_associated_df_sorted <- tf_hotspot_enh_associated_df[order(tf_hotspot_enh_associated_df$sites_percent, decreasing=TRUE)]
top_25_perc_for_enh_assoc <- as.numeric(quantile(tf_hotspot_enh_associated_df_sorted$sites_percent, 0.75)) ## or, use [[]] to extract the value

### Filetering the promoter associated region:
tf_hotspot_prom_associated <- tf_hotspot_region_based_file[tf_hotspot_region_based_file$state_anno == "prom_assoc"]
tf_hotspot_prom_associated$sites_percent <- (tf_hotspot_prom_associated$count/total_tf_hotspot_sites_prom)*100
tf_hotspot_prom_associated_sorted <- tf_hotspot_prom_associated[order(tf_hotspot_prom_associated$sites_percent, decreasing=TRUE)]
top_25_perc_for_prom_assoc <- as.numeric(quantile(tf_hotspot_prom_associated_sorted$sites_percent, 0.75))

### For whole genome - tf hotspot region based file:
output_dir <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/plots"
pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis.pdf"))

factor(tf_hotspot_region_based_file$tf_name) %>% levels
tf_hotspot_region_based_file$tf_name <- factor(tf_hotspot_region_based_file$tf_name, levels=tf_hotspot_file_sorted$tf_name)

ggplot(tf_hotspot_region_based_file, aes(x=tf_name, y=sites_percent)) + geom_bar(stat='identity') + 
          theme_gray() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="TF names", y="% of TF hotspot sites(n=12927) occupied")+
          theme(axis.text.y = element_text(size=2.0),
          plot.title=element_text(size=14, face="bold", hjust = 0.6)) +
          coord_flip()
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()

### For whole genome:
output_dir <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/plots"
pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis.pdf"))

### To maintain the order of TF in ggplots:
factor(tf_hotspot_file_sorted$tf_name) %>% levels
tf_hotspot_file_sorted$tf_name <- factor(tf_hotspot_file_sorted$tf_name, levels=tf_hotspot_file_sorted$tf_name)

ggplot(tf_hotspot_file_sorted, aes(x=tf_name, y=sites_percent)) + geom_bar(stat='identity') + 
          theme_gray() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="TF names", y="% of TF hotspot sites(n=12927) occupied")+
          theme(axis.text.y = element_text(size=2.0),
    			plot.title=element_text(size=14, face="bold", hjust = 0.6)) +
          coord_flip()
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()

pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis_top25_percent.pdf"))
ggplot(subset(tf_hotspot_file_sorted, sites_percent >= top_25_perc_whole_genome_assoc), aes(sites_percent, tf_name)) +
         geom_segment(aes(x = 0, y = tf_name, xend = sites_percent, yend = tf_name), color = "grey50") +
         geom_point() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(y="TF names", x="% of TF hotspot sites(n=12927) occupied")+
          theme(axis.text.y = element_text(size=5.0))
         #scale_x_continuous(labels=percent)
dev.off()

### For prom associated
output_dir = "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/plots"
pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis_prom_assoc.pdf"))

### To maintain the order of TF in ggplots:
factor(tf_hotspot_prom_associated_sorted$tf_name) %>% levels
tf_hotspot_prom_associated_sorted$tf_name <- factor(tf_hotspot_prom_associated_sorted$tf_name, levels=tf_hotspot_prom_associated_sorted$tf_name)

ggplot(tf_hotspot_prom_associated_sorted, aes(x=tf_name, y=sites_percent)) + geom_bar(stat='identity') + 
          theme_gray() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="TF names", y="% of Prom associated TF hotspot sites(n=6700) occupied")+
          theme(axis.text.y = element_text(size=2.0),
    			plot.title=element_text(size=14, face="bold", hjust = 0.6)) +
          coord_flip()
         # scale_y_continuous(labels=percent) + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()


pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis_top25_percent_prom_assoc.pdf"))
ggplot(subset(tf_hotspot_prom_associated_sorted, sites_percent >= top_25_perc_for_prom_assoc), aes(sites_percent, tf_name)) +
         geom_segment(aes(x = 0, y = tf_name, xend = sites_percent, yend = tf_name), color = "grey50") +
         geom_point() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(y="TF names", x="% of Prom associated TF hotspot sites(n=6700) occupied")+
          theme(axis.text.y = element_text(size=5.0))
         #scale_x_continuous(labels=percent)
dev.off()

### For enh associated
output_dir = "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/cobinding_analysis/plots"
pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis_enh_assoc.pdf"))

### To maintain the order of TF in ggplots:
factor(tf_hotspot_enh_associated_df_sorted$tf_name) %>% levels
tf_hotspot_enh_associated_df_sorted$tf_name <- factor(tf_hotspot_enh_associated_df_sorted$tf_name, levels=tf_hotspot_enh_associated_df_sorted$tf_name)

ggplot(tf_hotspot_enh_associated_df_sorted, aes(x=tf_name, y=sites_percent)) + geom_bar(stat='identity') + 
          theme_gray() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(x="TF names", y="% of Enh associated TF hotspot sites(n=6227) occupied")+
          theme(axis.text.y = element_text(size=2.0),
    			plot.title=element_text(size=14, face="bold", hjust = 0.6)) +
          coord_flip()
         # scale_y_continuous("Value") + 
         # scale_x_continuous("Distance from Center (bp)", limits=c(-20, 20))

dev.off()

pdf(file.path(output_dir,"TF_hotspot_cobinding_analysis_top25_percent_enh_assoc.pdf"))
ggplot(subset(tf_hotspot_enh_associated_df_sorted, sites_percent >= top_25_perc_for_prom_assoc), aes(sites_percent, tf_name)) +
         geom_segment(aes(x = 0, y = tf_name, xend = sites_percent, yend = tf_name), color = "grey50") +
         geom_point() +
          #facet_grid(variable~., scales = "free_y", labeller = as_labeller(meth_table_names) ) +
          labs(y="TF names", x="% of Enh associated TF hotspot sites(n=6227) occupied")+
          theme(axis.text.y = element_text(size=5.0))
         #scale_x_continuous(labels=percent)
dev.off()

