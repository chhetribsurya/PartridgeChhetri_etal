library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
library(tidyr)

output_dir  <- "/Users/suryachhetri/Dropbox/for_genemodels"
input_file <- file.path(output_dir, "ideas_stacked_prom_enh_piestate_dist_barplot_data_test.txt")

# Read input file:
input_df <- fread(input_file, sep="\t")
input_df$log2peak_count <- log2(input_df$total_count)
enh_assoc <- input_df[input_df$ideas_state_new == "Enhancer associated"]
prom_assoc <- input_df[input_df$ideas_state_new == "Promoter associated"]

prom_assoc <- prom_assoc %>% select(tf_name, ideas_state_new, bind_fraction, log2peak_count)
enh_assoc <- enh_assoc %>% select(tf_name, ideas_state_new, bind_fraction, log2peak_count)
merged_df <- merge(prom_assoc, enh_assoc, by.x=c("tf_name","log2peak_count"), by.y=c("tf_name", "log2peak_count"))
merged_df$tf_name <- merged_df$tf_name %>% str_replace("_human|\\[FLAG\\]", "")

# write.table(merged_df, file.path(output_dir, "revision_2", "prom_enh_fraction_binding_scatterplot_data.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

pdf(file.path(output_dir, "revision_2", "prom_enh_fraction_binding_scatterplot_1.pdf"))
plot <- ggplot(merged_df, aes(x=bind_fraction.y, y=bind_fraction.x, color=log2peak_count)) + 
        geom_point() + theme_bw() +
        xlab("Fraction in Enhancer assoc regions") + ylab("Fraction in Promoter assoc regions") + theme_bw()+ 
        ggtitle("Binding Fraction of Factors across Cis-Regulatory Regions") + 
        #guides(fill=guide_legend(title="Log2(peak_count)")) +
        theme(
        axis.text.x = element_text(size=9, face="bold"),
        axis.text.y = element_text(size=9, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust = 0.6)
        )

finalplot <- plot + geom_text(data = subset(merged_df, bind_fraction.x > 0.7 |  bind_fraction.y > 0.7 ), aes(x = bind_fraction.y, y =bind_fraction.x, label = tf_name), 
    hjust= -0.15, size = 1, color="black") #check_overlap=TRUE)

print(finalplot)
dev.off()

# #coord_flip() +
# ggtitle('Binding fraction of 208 factors across Cis-regulatory regions') +
# scale_x_discrete(name="Transcription Factors") +  scale_y_continuous(name="Binding Fraction") +
# #scale_fill_manual(name="IDEAS Anno",values=["red","orange"]) +
# guides(fill=guide_legend(title="Annotation")) +
# geom_hline(aes(yintercept=0.5), size=0.55, color="black", linetype="dashed") + 
# theme_minimal() + 
# theme(
# axis_text_y = element_text(size=1.3), #theme(axis_text_x=element_text(angle=45))
# plot_title = element_text(size=9, face="bold", hjust = 0.6),
# legend_title = element_text(size=8, face="bold") 
# ) +
# facet_wrap('~ideas_state_new')



# barplot <- barplot + geom_text(data = subset(meth_tf_df, meth_percent > 0.15), aes(x = celltype, y =meth_percent, label = tf_name), 
#     hjust= -0.15, size = 2, check_overlap = TRUE)


