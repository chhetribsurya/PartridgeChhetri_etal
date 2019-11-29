library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(gtools)
library(gplots)


### Barplot without TPM for TF count in each sp_gene bin:
# weights <- ifelse(pcaOutput2$PC1 < -5 & abs(pcaOutput2$PC2) > 10, 2, 1)

input_file <- "~/Dropbox/encode_3/hepg2_liver_gene_sp_analysis/files_2kb_from_tss/final_tf_binding_with_gene_bins_merged_barplot_data.bed"
tf_count_sp_gene_df <- fread(input_file, sep="\t", header= TRUE)

tf_count_sp_gene_df <- tf_count_sp_gene_df %>% 
                        data.frame %>%
                        arrange(tf_counts)

write_file = "~/Dropbox/Liver_hepg2_promoters_barplot_with_TPM_scale_data.txt"
write.table(tf_count_sp_gene_df, write_file, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE )

output_file_name <- paste0("~/Dropbox", "/", "final_tf_count_with_sp_gene_bin_barplot_test.pdf")                
pdf(output_file_name)
# test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))
reordered_gene_name = factor(tf_count_sp_gene_df$gene_name, levels=tf_count_sp_gene_df$gene_name)

barplot <- ggplot(tf_count_sp_gene_df, aes(x=reordered_gene_name, y=tf_counts)) + 
      geom_bar(stat="identity", fill="grey") +
    #geom_boxplot(aes(fill=tf_category),outlier.shape = NA) +
    ylab("Transcription Factor Counts") + xlab("Liver-HepG2 specific and highly expressed genes") + 
    theme_bw() + 
    ggtitle("Total TF counts across each genes") + 
    theme(
    axis.text.y = element_text(size=6, face="bold", color="black" ),
    plot.title=element_text(size=14, face="bold", hjust = 0.6),
    legend.title=element_text(face="bold")
    ) +
    #geom_text(aes(label=tf_counts), color="black", size=2.5, hjust = -0.1 ) +
    theme(axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust = 1)) +
    coord_flip() + 
    #scale_fill_gradient(low = "lightblue", high= "steelblue")+
    #scale_fill_manual(values=custom_col )+
    guides(fill=guide_legend(title="HepG2 TPM scale")) + 
    scale_y_continuous() 

print(barplot)
