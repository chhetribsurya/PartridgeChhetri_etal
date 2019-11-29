
library(data.table)
library(dplyr)
library(ggplot2)

# CEBPA
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "CEBPA.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "CEBPA.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "CEBPA_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("CEBPA ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# CEBPG
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "CEBPG.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "CEBPG.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "CEBPG_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("CEBPG ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# FOXA1
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "FOXA1.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "FOXA1.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "FOXA1_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("FOXA1 ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# FOXA2
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "FOXA2.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "FOXA2.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "FOXA2_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("FOXA2 ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# FOXA3
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "FOXA3.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "FOXA3.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "FOXA3_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("FOXA3 ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()

# HNF4A
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "HNF4A.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "HNF4A.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "HNF4A_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("HNF4A ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()

# HNF4A
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "HNF4A.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "HNF4A.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "HNF4A_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("HNF4A ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# HNF4G
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "HNF4G.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "HNF4G.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "HNF4G_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("HNF4G ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# KLF16
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "KLF16.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "KLF16.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "KLF16_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("KLF16 ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# RARA
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "RARA.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "RARA.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "RARA_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("RARA ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()


# SOX13_iso1
output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/barbara_hotnonhot_peaks_ChIPsignal"
read_df1 <- fread(file.path(output_dir, "SOX13_iso1.high.id.peaks"), sep="\t")
read_df2 <- fread(file.path(output_dir, "SOX13_iso1.low.id.peaks"), sep="\t")

# Column 7 shows ChIP signal:
ks_test_hotnonhot = ks.test(read_df1$V7, read_df2$V7) # p-value = 5.966e-11

pdf(file.path(output_dir, "SOX13_iso1_hot_nonhot_ChIPsignal_boxplot.pdf"))
boxplot(read_df1$V7, at=1, col=c("red2"), xlim=c(0, 3), main=paste0("SOX13_iso1 ChIP signal Distribution", "(KS-test < 2.2e-16)"))
boxplot(read_df2$V7, at=2, col=c("steelblue"), add=TRUE)
legend("topright", legend = c("HOTsites", "NonHOTsites"), col=c("red2", "steelblue"), pch = 15, bty = "n", pt.cex = 2, cex = 0.9,  horiz = F)
axis(1, at=seq(1,2,1), labels = c("HOT Sites(>70)", "NonHOT Sites(2-10)") , tick=FALSE , cex=0.3)
dev.off()











