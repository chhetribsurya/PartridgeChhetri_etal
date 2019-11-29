library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)
#library(xlsx)


args <-  commandArgs(TRUE)
input_file <- args[1]
tf_name <- args[2]
out_dir_name <- args[3]

#avg_df <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/wgbs_tf_intersect/final_wgbs_tf_intersect_average_meth.bed", sep="\t", header=TRUE)
xls_df <- read.xlsx2("~/Dropbox/Encode_full_hepg2_datasets_DBF_CR.xls", sheetIndex=4)
avg_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_intersect_average_meth.bed", sep="\t", header=TRUE)
meth_tf_df <- fread("~/Dropbox/encode_3/wgbs_tf_intersect_total/files/final_wgbs_tf_intersect.bed", sep="\t", header=TRUE)

names(avg_df) <- c("TF_name", "meth_percent")
names(meth_tf_df) <- c("TF_name", "index", "peak_chrom", "peak_start", "peak_end", "meth_percent")

target_vector <- avg_df$TF_name
xls_df_ordered <- xls_df[match(target_vector, xls_df$Target),]


# check if the order b/w 2 columns are equal: 
if(all(avg_df$TF_name == xls_df_ordered$Target))
{
   print("Column A(tf_name) and B(target) are identical")
}

xls_df_ordered <- xls_df_ordered %>% select(Target, Lab, Category)
avg_df$tf_category <-  xls_df_ordered$Category

cr_df <- avg_df %>% as.data.frame %>% filter(tf_category == "CR/CF")
dbf_df <- avg_df %>% as.data.frame %>% filter(tf_category == "DBF")
all_tf_df <- avg_df %>% as.data.frame

### or, do the following matching by ordering, but not necessary for this plot
### Not used in this plot, so could be removed if neeeded.
df1 <- xls_df %>%
    select(Target, Lab, Category) %>%
    arrange(Target) %>% as.data.table()

df2 <- avg_df %>% 
    arrange(TF_name) %>% as.data.table()

all.equal(df1,df2) # identical(df1,df2)

# Enrichment pvalue 0.001
read_file <- fread("~/Dropbox/encode_3/unique_hepg2_analysis_total/files/final_tf_enrichment_pval_corr_heatmap_data_0.001_state_filtered.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_enrichment_heatmap_0.001.png")				
output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_enrichment_heatmap_0.001pval.pdf")        
png(output_file_name, width = 7, height = 7, units = 'in', res = 600)
pdf(output_file_name)


summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(10) )
par(cex.main=0.7)
heatmap.2(mat_data,
  #dendrogram = "none",
  # Rowv = TRUE, 
  # Colv = TRUE,
  main = "HepG2 specific TF enrichment(<0.001)", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=improved_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  na.color="grey",
  breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

dev.off()



# Enrichment pvalue 0.01
read_file <- fread("~/Dropbox/encode_3/unique_hepg2_analysis_total/files/final_tf_enrichment_pval_corr_heatmap_data_0.01_state_filtered.bed", sep="\t", header=TRUE)
read_df <- as.data.frame(read_file)
read_df[is.na(read_df)] <- 0
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_enrichment_heatmap_0.01.png")        
output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_enrichment_heatmap_0.01_pval.pdf")        
png(output_file_name, width = 7, height = 7, units = 'in', res = 600)
pdf(output_file_name)

summary(read_df) #max_val = 325
improved_col <- c("#808080",greenred(10) )
par(cex.main=0.6)
heatmap.2(mat_data,
  #dendrogram = "none",
  #Rowv = TRUE, 
  #Colv = TRUE,
  main = "HepG2 specific TF enrichment(<0.01)", 
  xlab = "IDEAS States",
  ylab = "Transcription Factors",
  #col=greenred(10), 
  col=improved_col, 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.15,
  cexCol = 0.6,
  cex.main=0.3,
  na.color="grey",
  breaks=c(c(0,0.9),seq(1,325,(325-1)/9))
  #scale="column", tracecol="#303030"
  )    

dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count_heatmap.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "Tf_cobinding_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "TF Cobind Distribution", 
  xlab = "TF cobind",
  ylab = "TF name",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()

#### Unique hepg2 ideas tf analysis:
#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "TF name",
  ylab = "IDEAS state",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap_qval.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap_qval.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "IDEAS state",
  ylab = "TF name enriched with qvalue: < .001",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    

dev.off()


#### tf-tf cobinding analyis:
library(gplots)
#read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/tf_cobinding_analyis/final_cobinding_peak_count.bed", sep="\t", header=TRUE)
read_file <- fread("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/unique_hepg2_analysis/final_tf_ideas_peak_count_heatmap_2.bed", sep="\t", header=TRUE)
#heatmap_df <- as.matrix(read_df)
read_df <- as.data.frame(read_file)
rnames <- read_df[,1]
data <- read_df[,2:ncol(read_df)]
mat_data <- as.matrix(data)
rownames(mat_data) <- rnames
rownames(mat_data)
colnames(mat_data)


output_file_name <- paste0("~/Dropbox", "/", "hepg2_specific_tf_distribution_heatmap.pdf")				
pdf(output_file_name)

heatmap.2(mat_data,
  main = "HepG2 specific TF Distribution", 
  xlab = "TF name",
  ylab = "IDEAS state",
  col=greenred(10), 
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(6,8),      # widens margins around plot)
  cexRow = 0.35,
  cexCol = 0.35,
  na.color="grey"
  #scale="column", tracecol="#303030"
  )    
dev.off()




