#!/usr/bin/env python

library(data.table)
library(dplyr)

############################################

## Plot Redo : HOT motif permutation plot:
## Figure 6b:

#############################################

output_dir <- "/Users/suryachhetri/Dropbox/for_chris/hot_permutation_plot/redo1" 
output_dir <- "/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/HepG2full/cisbpmotif_id_extract" 

read_df <- fread(file.path(output_dir, "Hot_permutation_motif_percent.final.redo.srt.filt.txt"), sep="\t")

#calculate mean, min and max for each x-value
library(plyr)

read_df1 <- read_df[,c("tf_count", "one_motif_perc")]
df_final1 <- ddply(read_df1,.(tf_count),function(df) c(mean=mean(df$one_motif_perc),min=min(df$one_motif_perc),max=max(df$one_motif_perc), sdev=sd(df$one_motif_perc)))

read_df2 <- read_df[,c("tf_count", "two_motif_perc")]
df_final2 <- ddply(read_df2,.(tf_count),function(df) c(mean=mean(df$two_motif_perc),min=min(df$two_motif_perc),max=max(df$two_motif_perc), sdev=sd(df$two_motif_perc)))

read_df3 <- read_df[,c("tf_count", "three_motif_perc")]
df_final3 <- ddply(read_df3,.(tf_count),function(df) c(mean=mean(df$three_motif_perc),min=min(df$three_motif_perc),max=max(df$three_motif_perc), sdev=sd(df$three_motif_perc)))

read_df4 <- read_df[,c("tf_count", "true_hotsite_perc")]
df_final4 <- ddply(read_df4,.(tf_count),function(df) c(mean=mean(df$true_hotsite_perc),min=min(df$true_hotsite_perc),max=max(df$true_hotsite_perc), sdev=sd(df$true_hotsite_perc)))

read_df0 <- read_df[,c("tf_count", "zero_motif_perc")]
df_final0 <- ddply(read_df0,.(tf_count),function(df) c(mean=mean(df$zero_motif_perc),min=min(df$zero_motif_perc),max=max(df$zero_motif_perc), sdev=sd(df$zero_motif_perc)))


# 0 motif:
pdf(file.path(output_dir, "HOTSITE_permutation_plot.pdf"))

# 1 motif 
avg <- df_final1$mean
sdev <- df_final1$sdev
x <- df_final1$tf_count

plot(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    ylim=range(c(0, 100)),
    #xaxt = "n",
    pch=19, xlab="TF bound", ylab="Mean % of HOTSITES (+/- SD)",
    main="HOTSITE Permutation Plot", col = "red", xaxt="none", yaxt="none"

)

# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

# 2 motif:
avg <- df_final2$mean
sdev <- df_final2$sdev
x <- df_final2$tf_count
points(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    pch=19, col = "blue"
)
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

# 3 motif:
avg <- df_final3$mean
sdev <- df_final3$sdev
x <- df_final3$tf_count
points(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    pch=19, col = "green"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

# True hotsite perc
avg <- df_final4$mean
sdev <- df_final4$sdev
x <- df_final4$tf_count

points(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    #ylim=range(c(0, 100)),
    pch=19, col = "steelblue"
)

# # hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

#axis(1, labels=c(10,seq(10,140,10)), at=c(10,seq(10,140,10)))
axis(1, c(seq(20,160,20)), font=1)
axis(2, c(seq(0,100,20)), font=1)
legend("left", legend=c("% true HOT sites", "% with >=1 motif","% with >=2 motif", "% with >=3 motif"), pch=19, col=c("steelblue", "red", "blue", "green"))
# legend("left", legend=c("% with = 0 motif", "% with >=1 motif","% with >=2 motif", "% with >=3 motif"), pch=19, col=c("grey", "red", "blue", "green"))

dev.off()


########################################
# Plot without True HOT motif fraction:
##########################################

pdf(file.path(output_dir, "HOTSITE_permutation_plot_without_trueHotmotif.pdf"))

# 1 motif 
avg <- df_final1$mean
sdev <- df_final1$sdev
x <- df_final1$tf_count


plot(x, avg,
    ylim=range(c(avg-sdev, avg+sdev)),
    #ylim=range(c(0, 100)),
    #xaxt = "n",
    pch=19, xlab="TF bound", ylab="Mean % of HOTSITES (+/- SD)",
    main="HOTSITE Permutation Plot", col = "red", xaxt="none", yaxt="none"

)

# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

# 2 motif:
avg <- df_final2$mean
sdev <- df_final2$sdev
x <- df_final2$tf_count
points(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    pch=19, col = "blue"
)
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

# 3 motif:
avg <- df_final3$mean
sdev <- df_final3$sdev
x <- df_final3$tf_count
points(x, avg,
    #ylim=range(c(avg-sdev, avg+sdev)),
    pch=19, col = "green"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)

#axis(1, labels=c(10,seq(10,140,10)), at=c(10,seq(10,140,10)))
axis(1, c(seq(20,160,20)), font=1)
axis(2, c(seq(0,100,20)), font=1)
legend("topright", legend=c("% with >=1 motif","% with >=2 motif", "% with >=3 motif"), pch=19, col=c("red", "blue", "green"))
# legend("topright", legend=c("% with = 0 motif", "% with >=1 motif","% with >=2 motif", "% with >=3 motif"), pch=19, col=c("grey", "red", "blue", "green"))

dev.off()
