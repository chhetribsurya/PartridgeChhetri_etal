library(tidyverse)
library(readxl)
library(dplyr)
library(data.table)

xlsx_example <- "/Users/suryachhetri/Dropbox/for_chris/concordant_disconcordant/20190618_temp_supp_table_4.xls"
excel_sheets(xlsx_example)
cisbp_scan_tomfile_xls <- read_excel(xlsx_example, sheet="Supplementary_table_4") %>% as.data.frame
cisbp_scan_tomfile_new <- cisbp_scan_tomfile_xls[,c(2:7)] 
#names(cisbp_scan_tomfile_new) <- c("TF_NAME", "MOTIF_ID", "Closest_target", "anno_redo")

# Custom inhouse script (off peak distance):
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset_nearest.txt"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset_nearest_redo.txt"
tomtom_df = fread(file.path(output_dir, cisbp_scan_tomfile), sep="\t")
names(tomtom_df)[names(tomtom_df) == 'tom_E-value_y'] <- "tom_E_value_final"
names(tomtom_df)[names(tomtom_df) == 'cmo_E-value'] <- "cmo_E_value"
#names(tomtom_df) <- c("motif_fraction", "tom_Evalue")

output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_offpeak_dist_nearest.pdf")       
output_file_name <- file.path(output_dir, "Distant_motif_similarity_plot_with_293_motifs_offpeak_dist_nearest_redo.pdf")       
pdf(output_file_name)

# Excluding "no_match"
tomtom_df_subset = tomtom_df[!tomtom_df$anno == "no_match",]
match_subset = tomtom_df[tomtom_df$anno == "match",]
mismatch_subset = tomtom_df[tomtom_df$anno == "mismatch",]
nomatch_subset = tomtom_df[tomtom_df$anno == "no_match",]
ks_test = ks.test(match_subset$adj_offset_dist,mismatch_subset$adj_offset_dist)
ggplot(tomtom_df_subset, aes(y=as.factor(tf_motif_id), x=adj_offset_dist, label = tf_motif_id )) +
        geom_point(aes(color=anno_new), size=0.9) + geom_text(check_overlap = TRUE, size=1.5) +
        scale_color_manual(values=c("green","orange","red")) +
        #geom_rug(color="darkgray", sides="b") +
        geom_vline(xintercept=0,linetype=2) +
        #geom_vline(xintercept=14,linetype=2) +
        labs(x="Offpeak median Distance", y="DNA associated factors", color='') + # xlabel, ylabel and legend title
        #scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        #scale_y_continuous(breaks=seq(0,50,5)) +
        ggtitle(paste("Kolmogorov-Smirnov Test (pval) :", 2.2e-16)) + #ks_test$p.value
        theme_bw()+ theme(plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.y = element_text(size=0),
             legend.title=element_text(face="bold")
        ) + facet_wrap(~anno_new) 
        
#       guides(fill=guide_legend(title="Annotation")) + theme_bw()
            # ylim(0,1)
dev.off()

# Novel motif no match subset:
output_file_name <- file.path(output_dir, "Offpeak_dist_no_match_motifs.pdf")       
pdf(output_file_name)
ggplot(nomatch_subset, aes(y=as.factor(tf_motif_id), x=adj_offset_dist, label = tf_motif_id )) +
        geom_point(color="green", size=0.9) + geom_text(check_overlap = TRUE, size=2.5) +
        scale_color_manual(values=c("green","orange","red")) +
        #geom_rug(color="darkgray", sides="b") +
        geom_vline(xintercept=0,linetype=2) +
        #geom_vline(xintercept=14,linetype=2) +
        labs(x="Offpeak median Distance", y="DNA associated factors", color='') + # xlabel, ylabel and legend title
        #scale_x_continuous(breaks=seq(0.05,1,0.1)) +
        #scale_y_continuous(breaks=seq(0,50,5)) +
        ggtitle("Motif Location Distribution") + #ks_test$p.value
        theme_bw()+ theme(plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.y = element_text(size=0),
             legend.title=element_text(face="bold")
        ) + facet_wrap(~anno_new) 
        
#       guides(fill=guide_legend(title="Annotation")) + theme_bw()
            # ylim(0,1)
dev.off()


#############################################
# For stacked histogram plot with base R:
# Concordance and disconcordance calls based
# plots
#############################################

library(tidyverse)
library(readxl)
library(dplyr)
library(data.table)

### Newly re-annotated xls file read:
xlsx_example <- "/Users/suryachhetri/Dropbox/for_chris/concordant_disconcordant/20190823_temp_supp_table_4.xls"
# excel_sheets(xlsx_example)

cisbp_scan_tomfile_xls <- read_excel(xlsx_example, sheet="Supplementary_table_4") %>% as.data.frame
cisbp_scan_tomfile_new <- cisbp_scan_tomfile_xls[,c(2,3,4,7)] 
names(cisbp_scan_tomfile_new) <- c("TF_NAME", "MOTIF_ID", "Closest_target", "anno_redo")
cisbp_scan_tomfile_new <- na.omit(cisbp_scan_tomfile_new)

### Previous offset distance containing file:
output_dir <- "/Users/suryachhetri/Dropbox/for_chris/concordant_disconcordant"
input_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
cisbp_scan_tomfile <- "Novel_motifs_combined_tomtom_cisbp_scan_df_top50_peakoffset_nearest_redo.txt"
tomtom_df <- fread(file.path(input_dir, cisbp_scan_tomfile), sep="\t") %>% as.data.frame
tomtom_df <- tomtom_df[,c("TF_NAME", "MOTIF_ID", "adj_offset_dist", "anno_new")]
names(tomtom_df) <- c("TF_NAME", "MOTIF_ID", "adj_offset_dist", "anno_previous")

# names(tomtom_df)[names(tomtom_df) == 'tom_E-value_y'] <- "tom_E_value_final"
# names(tomtom_df)[names(tomtom_df) == 'cmo_E-value'] <- "cmo_E_value"
# tomtom_df %>% filter(anno_previous=="discordance") %>% nrow
# [1] 165
# tomtom_df %>% filter(anno_previous=="concordance") %>% nrow
# [1] 103
# tomtom_df %>% filter(anno_previous=="no_match") %>% nrow
# [1] 26

# Merge both dataframe xls based re-annotated and old annotation of motif:
merged_tomtom_df <- merge(cisbp_scan_tomfile_new, tomtom_df, by.x=c("TF_NAME", "MOTIF_ID"), by.y=c("TF_NAME", "MOTIF_ID"), all.x=TRUE)
merged_tomtom_df <- na.omit(merged_tomtom_df)
# merged_tomtom_df %>% filter(is.na(anno_previous))

##########################################
# or directly read the file pre-produced:
xlsx_example <- "/Users/suryachhetri/Dropbox/for_chris/concordant_disconcordant/Merged_df_motif_offset_dist.xlsx"
# excel_sheets(xlsx_example)

merged_tomtom_df <- read_excel(xlsx_example) %>% as.data.frame
merged_tomtom_df <- na.omit(merged_tomtom_df)
merged_tomtom_df$adj_offset_dist <- round(merged_tomtom_df$adj_offset_dist)

##########################################

#merged_tomtom_df_subset = merged_tomtom_df[!merged_tomtom_df$anno == "no_match",]
concordant_df <- merged_tomtom_df %>% filter(redo_anno=="concordant")
discordant_df <- merged_tomtom_df %>% filter(redo_anno=="discordant")
nomatch_df <- merged_tomtom_df %>% filter(redo_anno=="no_match")
ks_test = ks.test(concordant_df$adj_offset_dist,discordant_df$adj_offset_dist)

vec1 <- concordant_df$adj_offset_dist
vec2 <- discordant_df$adj_offset_dist
vec3 <- nomatch_df$adj_offset_dist
bin_width <- 3

output_file_name <- file.path(output_dir, "Median_offset_dist_plot_stacked.pdf")       
pdf(output_file_name)
hist(concordant_df$adj_offset_dist, col=rgb(0,0,1,0.25), xlim=c(0,60), ylim=c(0,25), breaks=seq(min(vec1),max(vec1)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
hist(discordant_df$adj_offset_dist, xlim=c(0,60), ylim=c(0,25), col=rgb(1,0,0,0.25),breaks=seq(min(vec2),max(vec2)+2, bin_width), add=T)
abline(v=c(median(vec1),median(vec2)), col = c(rgb(0,0,1,0.25), col=rgb(1,0,0,0.25)), lty=c(2,2), lwd=c(4, 4))
legend("topright", legend=c("Concordant","Discordant"), pch=15, col=c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)))
dev.off()

output_file_name <- file.path(output_dir, "Median_offset_dist_plot_concordant.pdf")       
pdf(output_file_name)
hist(concordant_df$adj_offset_dist, col=rgb(0,0,1,0.25), xlim=c(0,50), ylim=c(0,25), breaks=seq(min(vec1),max(vec1)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec1), lty=2, lwd=2)
dev.off()

output_file_name <- file.path(output_dir, "Median_offset_dist_plot_discordant.pdf")       
pdf(output_file_name)
hist(discordant_df$adj_offset_dist, col="orange", xlim=c(0,50), ylim=c(0,25), breaks=seq(min(vec2),max(vec2)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec2), lty=2, lwd=2)
dev.off()

output_file_name <- file.path(output_dir, "Median_offset_dist_plot_no_match.pdf")       
pdf(output_file_name)
hist(nomatch_df$adj_offset_dist, col="green", xlim=c(0,50), ylim=c(0,25), breaks=seq(min(vec3),max(vec3)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec3), lty=2, lwd=2)
dev.off()


# > mean(vec1)
# [1] 14.76068
# median(vec1)
# [1] 12
# median(vec2)
# [1] 23
# median(vec3)
# [1] 9

# > c(median(vec1),median(vec2), median(vec3))
# [1]  9 24 11
# mean(vec1)
# [1] 10.26699
# mean(vec2)
# [1] 28.07879
# mean(vec3)
# [1] 20.73077

# > median(vec1) # concordant
# [1] 12
# > median(vec2) # disconcordant
# [1] 23
# > median(vec3) # no_match
# [1] 9


###########################

# Figure 3b and supp 9 fig:

# Combining all 3 plots:

###########################


output_file_name <- file.path(output_dir, "Median_offset_dist_plot_all_combined_newupdated.pdf")       
pdf(output_file_name)

par(mfrow=c(2,2)) #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
hist(concordant_df$adj_offset_dist, col=rgb(0,0,1,0.25), xlim=c(0,max(concordant_df$adj_offset_dist)), ylim=c(0,25), breaks=seq(min(vec1),max(vec1)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec1), lty=2, lwd=2)

hist(discordant_df$adj_offset_dist, col="orange", xlim=c(0,100), ylim=c(0,25), breaks=seq(min(vec2),max(vec2)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec2), lty=2, lwd=2)

hist(nomatch_df$adj_offset_dist, col="green", xlim=c(0,100), ylim=c(0,25), breaks=seq(min(vec3),max(vec3)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
abline(v=median(vec3), lty=2, lwd=2)

hist(concordant_df$adj_offset_dist, col=rgb(0,0,1,0.25), xlim=c(0,100), ylim=c(0,25), breaks=seq(min(vec1),max(vec1)+2, bin_width), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
hist(discordant_df$adj_offset_dist, xlim=c(0,100), ylim=c(0,25), col="orange",breaks=seq(min(vec2),max(vec2)+2, bin_width), add=T)
hist(nomatch_df$adj_offset_dist, xlim=c(0,100), ylim=c(0,25), col="green", breaks=seq(min(vec2),max(vec2)+2, bin_width), add=T)
abline(v=c(median(vec1),median(vec2), median(vec3)), col = c(rgb(0,0,1,0.25), col="orange", "green"), lty=c(2,2,2), lwd=c(3, 3, 3))
legend("topright", legend=c("Concordant","Discordant", "no_match"), pch=15, col=c(rgb(0,0,1,0.25),"orange", "green"))
#abline(v=c(median(vec1),median(vec2), median(vec3)), col = c(rgb(0,0,1,0.25), col=rgb(1,0,0,0.25), "green"), lty=c(2,2,2), lwd=c(3, 3, 3))
#legend("topright", legend=c("Concordant","Discordant", "no_match"), pch=15, col=c(rgb(0,0,1,0.25),rgb(1,0,0,0.25), "green"))

dev.off()

