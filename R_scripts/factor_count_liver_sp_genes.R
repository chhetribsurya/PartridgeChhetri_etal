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
dev.off()


# Barplot with TPM value i.e with fill gradient:
output_file_name <- paste0("~/Dropbox", "/", "Liver_hepg2_promoters_barplot_with_TPM_scale.pdf")        
pdf(output_file_name)
# test <- data.frame(X = c("DBF", "CR/CF" ), Z = c(0.15, 0.15))
reordered_gene_name = factor(tf_count_sp_gene_df$gene_name, levels=tf_count_sp_gene_df$gene_name)

barplot <- ggplot(tf_count_sp_gene_df, aes(x=reordered_gene_name, y=tf_counts, fill=HepG2_RNAseq)) + 
    geom_bar(stat="identity") +
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
    scale_fill_gradient(low="lightblue", high= "darkblue")+
    #scale_fill_gradient2(low="grey", mid = "lightblue", high= "steelblue")+
    #scale_fill_manual(values=custom_col )+
    guides(fill=guide_legend(title="HepG2 TPM scale")) + 
    scale_y_continuous() 

print(barplot)
dev.off()

# Heatmap:
input_file <- "~/Dropbox/encode_3/hepg2_liver_gene_sp_analysis/files_2kb_from_tss/final_tf_binding_with_gene_bins_heatmap_data.bed"
#input_file <- "~/Dropbox/encode_3/tf_factor_count_with_sp_gene_bins/files/final_tf_sp_gene_piechart_combined_heatmap_data.bed"
read_file <- fread(input_file, sep="\t", header=TRUE)

read_df <- as.data.frame(read_file)
#read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "final_tf_count_with_sp_gene_bin_heatmap.pdf")				
pdf(output_file_name)

#green_col <-  colorRampPalette("green")
#red_col1 <-  colorRampPalette("red")

grey_col <-  colorRampPalette("grey")
red_col <- colorRampPalette(c("indianred1"))

#custom_col <- c(green_col(1), red_col1(1))
custom_col <- c(grey_col(1), red_col(1))

summary(read_df) #max_val = 325
#improved_col <- c("#808080",greenred(10) )
par(cex.main=0.55)
test <- heatmap.2(mat_data,
  #dendrogram = "none",
  # Rowv = TRUE, 
  # Colv = TRUE,
  #scale="none",
  main = "TF occupancy across 57 Liver-HepG2 specific/expressed genes", 
  xlab = "Gene names",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=custom_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.45,
  na.color="grey"
  #symm=F,symkey=T,symbreaks=T,
  #breaks=c(-1,0.8,1.2,3)

  #breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    
test
dev.off()


# Plot the strip heatmap of RPKM for genes in 
# the same order as the main/central heatmap:
order_of_heatmap <- mat_data[rev(test$rowInd), test$colInd]
gene_names <- colnames(order_of_heatmap)
tf_names <- rownames(order_of_heatmap)

heatmap_ordered_df <- tf_count_sp_gene_df[match(gene_names, tf_count_sp_gene_df$gene_name),]
melted_stacked_df <- melt(heatmap_ordered_df, id.vars="gene_name", measure.vars =c("HepG2_RNAseq", "Liver_RNAseq") )
melted_stacked_df$re_gene_name <- factor(melted_stacked_df$gene_name, levels=melted_stacked_df$gene_name)
casted_unstacked_df <- cast(melted_stacked_df, variable~re_gene_name) # but, not used for these analysis

output_file_name <- paste0("~/Dropbox", "/", "Liver_hepg2_promoters_heatmap_TPM_chart_1.pdf")       
pdf(output_file_name)

p <- ggplot(melted_stacked_df, aes(re_gene_name, variable)) +
  geom_tile(aes(fill = value), colour = "black") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  #scale_fill_gradient2(low=muted("blue"), high=muted("red"))+
  #scale_fill_gradient(low="red", high="yellow") +
  #scale_fill_gradient(low = "white", high = "red") + 
  labs(x="Genes", y="", title="") +
  geom_text(aes(label = round(value, 1)), size=0.6) +
  theme(axis.text.x = element_text(size = 2.5, angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(size = 4, colour = "black")) +
  #coord_flip() 
  coord_fixed(ratio = 1.5)
  #theme_dendro()
  #facet_wrap(~variable)
p
dev.off()


