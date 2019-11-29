library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(gtools)
library(gplots)


###############################################
#
# All foxa3, gatad2a, random sites gene exp:
#
###############################################

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/revision_2"
gene_exp_df <- fread(file.path(output_dir, "gatad2a_foxa3_nearest_tss_gene_expr.txt"), sep="\t", header=TRUE)
gene_exp_df$log2TPM <- log2(gene_exp_df$TPM+1)
#names(meth_tf_df) <- c("celltype", "tf_name", "TPM")

## Output file:
output_file_name <- file.path(output_dir, "gatad2a_foxa3_nearest_tss_gene_expr.pdf")             
pdf(output_file_name)


barplot <- ggplot(gene_exp_df, aes(x=annotation, y=log2(TPM) +1, fill=annotation)) + 
    geom_boxplot(width=0.25, outlier.colour = "red") + # outlier.shape = NA to remove the outliers 
    #geom_dotplot(data=gene_exp_df, binaxis='y', stackdir='centerwhole', dotsize=0.2, fill="black") +
    #stat_boxplot( aes(annotation, log2(TPM) +1), 
    #geom='errorbar', linetype=2, width=0.15, color="blue") +  
    xlab("Annotation") + ylab("HepG2 RNAseq Log2(TPM +1)") + theme_bw()+ 
    ggtitle("Gene Expression Distribution") + 
    #scale_y_continuous(limits=c(0,1)) +
    theme(
    axis.text.x = element_text(size=10, face="bold"),
    plot.title=element_text(size=14, face="bold", hjust = 0.6)
    )
#barplot <- barplot + geom_text(data = subset(meth_tf_df, log2(TPM) +1 > 0.15), aes(x = annotation, y =TPM, label = tf_name), 
#    hjust= -0.15, size = 2, check_overlap = TRUE)

print(barplot)
dev.off()

ks.test(x,y)



# Test Experimental boxplots:
gene_exp_df <- fread(file.path(output_dir, "gatad2a_foxa3_random_combined_gene_exp.txt"), sep="\t", header=TRUE)
gene_exp_df$log2TPM <- log2(gene_exp_df$TPM+1)

output_file_name <- file.path(output_dir, "gatad2a_foxa3_random_nearest_geneExp_boxplot.pdf")       
pdf(output_file_name)

color <- c("green", "red", "blue", "grey", "grey") # order the color based on categories
a <- boxplot(gene_exp_df, 
            col=color, 
            boxwex=0.55, 
            cex.axis=0.7, 
            ylab="HepG2 RNA-Seq (TPM)",
            xlab="Category",
            main="Binding Sites Nearest Gene Exp",
            font=2,
            names=c("GATAD2A_FOXA3", "FOXA3_only", "GATAD2A_only", "FOXA3_random", "GATAD2A_random"), 
            frame=TRUE, 
            outline=FALSE
            ) # outline=TRUE for outliers, #boxwex=0.4, main="title"
# plot(a)
text(cex=0.7, a$stats[nrow(a$stats) , ]+2 , paste("Median = ",c("15.18", "12.31", "13.93", "1.73", "1.94"), sep="")  )
#mtext("KS_Test Pval:  GATAD2A_FOXA3 vs FOXA3_only = 5.22*10-19
#                      GATAD2A_FOXA3 vs GATAD2A_only = 1.52*10-06
#                      GATAD2A_FOXA3 vs FOXA3_random = 0.00 ", cex=0.4, side=3)
# axis(1, labels = c("gatad2a_foxa3", "foxa3_only", "gatad2a", "foxa3", "gatad2a"), at=c(1:5),las=2,cex.axis=2)

legend("topright", legend = " 
                      GATAD2A_FOXA3 vs FOXA3_only = 5.22 * 10-19
                      GATAD2A_FOXA3 vs GATAD2A_only = 1.52 * 10-06
                      GATAD2A_FOXA3 vs FOXA3_random = 0.00 
                      GATAD2A_FOXA3 vs GATAD2A_random = 0.00 ", 
                      cex=0.45, bty = "n", col=c("red"),
                      text.font=2,
                      text.col="red",
                      title="KS_Test Pvalue",
                      title.col="black") 

                  # col = c("red" , "green", "blue") , bty = "n", pch=20 , pt.cex = 3, cex = 1, horiz = FALSE, inset = c(0.03, 0.1))
# boxplot(mat_data2 , col=ifelse(colnames(mat_data2)=="FOXA3" , "red",  "lightblue"), xlab="no box", ylim=c(0,max(mat_data2)), frame=F) #boxwex=0.4 , main=""
# boxplot(mat_data2 , col=ifelse(levels(mat_data2$tf_name)=="FOXA3" , rgb(0.1,0.1,0.7,0.5),  "lightblue"), xlab="no box", ylim=c(0,max(mat_data2)), frame=F) #boxwex=0.4 , main=""

dev.off()

