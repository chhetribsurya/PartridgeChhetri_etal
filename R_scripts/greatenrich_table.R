### Great analysis for the TF bound and unbound promoter nearby genes:
library(dplyr)
library(ggplot2)
library(ggExtra)
require(gridExtra)

input_dir <- "~/Dropbox/encode_3/tf_hotspots_prom_enh_categorised/files_d_100/merged_d_100/great_analysis/downloaded_final_data"
input_file <- file.path(input_dir, "All_TF_bound_unbound_promoter_nearby_genes_w_tpm_vals.txt")
output_dir <- file.path(input_dir, "plots")

## Provide the dir name(i.e sub dir) that you want to create under main dir:
#ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
if (!dir.exists(output_dir)){
	dir.create(output_dir)
} else {
	print("Dir already exists!")
}

tf_bound_unbound_final <- fread(input_file, sep="\t", header=TRUE)
tf_bound_unbound_final$log10_TPM = log10(tf_bound_unbound_final$Value)
tf_bound_unbound_final$log10_TPM[tf_bound_unbound_final$log10_TPM == -Inf] <- 0

tf_bound_unbound_final$log2_TPM = log2(tf_bound_unbound_final$Value)
tf_bound_unbound_final$log2_TPM[tf_bound_unbound_final$log2_TPM == -Inf] <- 0

output_file_name <- file.path(output_dir, "TF_bound_unbound_promoter_nearby_genes_tpm_analysis_boxplot.pdf")			
pdf(output_file_name)

barplot <-ggplot(tf_bound_unbound_final, aes(x=anno, y=log2_TPM, fill=anno)) + 
    geom_boxplot(width=0.4, outlier.colour = NA) + # outlier.shape = NA to remove the outliers 
    stat_boxplot( aes(anno, log2_TPM), 
    geom='errorbar', linetype=2, width=0.15, color="blue") + 
    annotate(geom="text", x=2.2, y=13.5, label = 'atop(bold("KS 2_Sample_Test"),"pvalue = 2.23e-150")',
            colour = "red", parse = TRUE) +
    xlab("Promoter associated sites") + ylab("Gene expression log10(TPM)") + theme_bw()+ 
    ggtitle("Gene expression distribution on promoter associated TF bound and unbound sites") + 
    scale_y_continuous() +
    guides(fill=guide_legend(title="Annotation")) +
    theme(
    axis.text.x = element_text(size=10),
    plot.title=element_text(size=10,hjust = 0.6, face="bold")
    ) 

print(barplot)

dev.off()
    
