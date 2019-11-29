
## For stacked histogram plot with base R:
## Local machine plotting:

library(data.table)
library(dplyr)
library(ggplot2)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment"
hot_read_df <- fread(file.path(output_dir, "Hotmotif_sites_histogram_data.txt"), sep="\t")
nonhot_read_df <- fread(file.path(output_dir, "NonHotmotif_sites_histogram_data.txt"), sep="\t")
random_read_df <- fread(file.path(output_dir, "Random_sites_histogram_data.txt"), sep="\t")

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm_peaklevel"
hot_read_df <- fread(file.path(output_dir, "Hotmotif_sites_histogram_data_tfbased_10.txt"), sep="\t")
nonhot_read_df <- fread(file.path(output_dir, "NonHotmotif_sites_histogram_data_tfbased_10.txt"), sep="\t")
random_read_df <- fread(file.path(output_dir, "Random_sites_histogram_data_tfbased_10.txt"), sep="\t")

df_10 <- hot_read_df %>% filter(percent_classifier == "top10_percent")
df_75 <- hot_read_df %>% filter(percent_classifier == "bottom75_percent")
df_10$num_bound_tfs %>% mean #2.888125
df_75$num_bound_tfs %>% mean

df_10 <- nonhot_read_df %>% filter(percent_classifier == "top10_percent")
df_75 <- nonhot_read_df %>% filter(percent_classifier == "bottom75_percent")
df_10$num_bound_tfs %>% mean #0.5306554
df_75$num_bound_tfs %>% mean

df_10 <- random_read_df %>% filter(percent_classifier == "top10_percent")
df_75 <- random_read_df %>% filter(percent_classifier == "bottom75_percent")
df_10$num_bound_tfs %>% mean #0.4252995
df_75$num_bound_tfs %>% mean

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/enhancer_random_nonhotsite_hotsite_factors_enrichment"
pdf(file.path(output_dir, "Hotmotif_sites_histogram_svm_clf_val_finplot_10.pdf"))
plot <- ggplot(hot_read_df, aes(x=num_bound_tfs, fill=percent_classifier))+ 
        geom_histogram(binwidth=1, position="identity", alpha=0.7) + 
        #ggtitle("Ranked Classifier-Weights Distribution") +
        theme_bw() +
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Number of Hotsites(>70 TFs cobound)") + xlab("Number of bound TFs with SVM classifier values(each site)") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="Ranked Classifier Value")) +
        theme(
        axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot.title = element_text(size=14, face="bold"),
        legend.title = element_text(size=8, face="bold"),
        legend.position = "top"
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
plot
dev.off()


pdf(file.path(output_dir, "NonHotmotif_sites_histogram_svm_clf_val_finplot_10.pdf"))
plot <- ggplot(nonhot_read_df, aes(x=num_bound_tfs, fill=percent_classifier))+ 
        geom_histogram(binwidth=1, position="identity", alpha=0.7) + 
        #ggtitle("Ranked Classifier-Weights Distribution") +
        theme_bw() + 
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Number of Hotsites (2-10 TFs cobound)") + xlab("Number of bound TFs with SVM classifier values(each site)") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="Ranked Classifier Value")) +
        theme(
        axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot.title = element_text(size=14, face="bold"),
        legend.title = element_text(size=8, face="bold"),
        legend.position = "top"
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
plot
dev.off()


pdf(file.path(output_dir, "Random_Genome_sites_histogram_svm_clf_val_finplot_10.pdf"))
plot <- ggplot(random_read_df, aes(x=num_bound_tfs, fill=percent_classifier))+ 
        geom_histogram(binwidth=1, position="identity", alpha=0.7) + 
        #ggtitle("Ranked Classifier-Weights Distribution") +
        theme_bw() + 
        # scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
        ylab("Number of Hotsites(>70 TFs cobound)") + xlab("Number of bound TFs with SVM classifier values(each site)") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="Ranked Classifier Value")) +
        theme(
        axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot.title = element_text(size=14, face="bold"),
        legend.title = element_text(size=8, face="bold"),
        legend.position = "top"
        # axis_text_y = element_text(size=1.3),
        # axis_text_y = element_text(size=1.3) 
        )
plot
dev.off()


library(plyr)
mu <- ddply(read_df, "percent_classifier", summarise, grp.mean=mean(num_bound_tfs))
head(mu)
plot + geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")

# #####################
# top5_df <- read_df %>% filter(percent_classifier=="top5_percent")
# bottom75_df <- read_df %>% filter(percent_classifier=="bottom75_percent")
# vec1 <- top5_df$num_bound_tfs
# vec2 <- bottom75_df$num_bound_tfs
# bin_width <- 5

# output_file_name <- file.path(output_dir, "Median_offset_dist_plot_stacked.pdf")       
# pdf(output_file_name)
# ggplot(read_df) + geom_hist
# hist(vec1, col=rgb(0,0,1,0.25), main="Offset Distribution", ylab="Number of Motifs", xlab="Median offset distance(bp)")
# hist(vec2, col=rgb(1,0,0,0.25), add=T)
# abline(v=c(median(vec1),median(vec2)), col = c(rgb(0,0,1,0.25), col=rgb(1,0,0,0.25)), lty=c(2,2), lwd=c(4, 4))
# legend("topright", legend=c("Concordant","Discordant"), pch=15, col=c(rgb(0,0,1,0.25),rgb(1,0,0,0.25)))
# dev.off()

# dev.off()

###########################################################################
###########################################################################


