
# R ggplot locally
library(data.table)
library(ggplot2)
library(dplyr)

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
plot_df = fread(file.path(output_dir, "Hotmotif_sites_problem2_histogram_data.txt"), sep="\t")
plot_df = fread(file.path(output_dir, "Hotmotif_sites_problem2_histogram_data_enhbound.txt"), sep="\t")
plot_df = fread(file.path(output_dir, "Hotmotif_sites_problem2_histogram_data_enhbound_random5000.txt"), sep="\t")
plot_df = fread(file.path(output_dir, "Hotmotif_sites_problem2_histogram_data_enhbound_all_24679.txt"), sep="\t")
out_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm"
clf_mean_count <- plot_df %>% group_by(percent_classifier) %>% summarise(mean_count=mean(num_bound_tfs)) %>% select(mean_count)

pdf(file.path(out_dir, "Hotmotif_sites_problem2_histogram_svm_clf_value_figure.pdf"))
pdf(file.path(out_dir, "Hotmotif_sites_problem2_histogram_svm_clf_value_figure_enh_bound_6_10_random2500.pdf"))
pdf(file.path(out_dir, "Hotmotif_sites_problem2_histogram_svm_clf_value_figure_enh_bound_6_10_random5000.pdf"))
pdf(file.path(out_dir, "Hotmotif_sites_problem2_histogram_data_enhbound_6_10_all_24679.pdf"))
plot <- ggplot(plot_df, aes(x=num_bound_tfs, fill=percent_classifier)) + 
        geom_histogram(binwidth=0.5) +
        geom_vline(linetype="dashed", color="red", xintercept=clf_mean_count$mean_count[1]) +
        geom_vline(linetype="dashed", color="steelblue", xintercept=clf_mean_count$mean_count[2]) +

        ggtitle("Ranked Classifier-Weights Distribution") +
        # theme_bw() +
        scale_x_continuous(breaks=seq(0,10)) +
        ylab("Number of Hotsites(6-10 TFs cobound)") + xlab("Number of bound TFs with SVM classifier values (each site)") +
        #scale_fill_manual(name="IDEAS Anno",values=["blueviolet","orange","green","red", "burlywood"]) +
        guides(fill=guide_legend(title="Ranked Classifier Value")) +
        theme(
        axis.title.y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        axis.title.x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
        plot.title = element_text(size=14, face="bold"),
        legend.title = element_text(size=8, face="bold"))
plot

dev.off()
