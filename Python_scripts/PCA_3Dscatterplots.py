import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import matplotlib; import matplotlib.pyplot as plt
import os; from os.path import join, expanduser
#from ggplot import *
from plotnine import *

##########################################################
# PCA analysis based on TFs binding at chromatin states
##########################################################

output_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
input_dir = "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"

#ideas_table <- read_excel(file.path(output_dir, "IDEAS_reannotation_table.xlsx" ))
ideas_table = pd.read_excel(join(input_dir, "IDEAS_reannotation_table.xlsx"))
ideas_table["final_ano"] = ideas_table["annotation"].str.replace(r"Cluster .* \(", "").str.replace(r"regions\)", "").str.replace(r"\)", "").str.replace(" ", "_").str.replace("__", "").str.replace("-", "_").str.replace("Body_", "Body")
ideas_table["label"] = ideas_table["final_ano"]
ideas_table["label"] = ideas_table["label"].replace({"Promoter_like":"red", "Enhancer_like":"orange", "Promoter_and_Enhancer_like":"#619CFF", "Heterochromatin_and_Repressed":"gray67", "CTCF_Bound":"purple3", "Active_Gene_Body":"green"})
ideas_table.to_csv("~/Dropbox/TFs_reAnnotation_file.txt", sep="\t", header=True, index=False)


input_file = expanduser("~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt")
df_orig = pd.read_csv(input_file, sep="\t") # original value containing file
df_new = pd.read_csv("~/Dropbox/TFs_reAnnotation_file.txt", sep="\t") # TF annotation file
df_merged = pd.merge(df_orig, df_new, left_on="tf_name", right_on="Target")

features = df_orig.columns # 
target = "annotation"
# target = "Category"

# Separating out the features and convert to numpy arrays:
x = df_merged.loc[:, features ].set_index("tf_name").values # tf_name is string so set_index()

# Separating out the target or labels: 
y = df_merged.loc[:, [target] ].values

# Standardizing the features:
# x = StandardScaler().fit_transform(x)
# Scaling the values:
# x = scale(x)

# PCA analysis (PCA projection to 2D/3D):
pca = PCA() # PCA(n_components=5)
principal_comp = pca.fit_transform(x) # X already normalized with fraction
print("original shape:   ", x.shape)
print("transformed shape:", principal_comp.shape) # eigen decomposed data(reduced to n dim)

# Transform to df:
principal_comp_df = pd.DataFrame(principal_comp)

# The anmount of variance that each PC explains
comp = pca.components_ # accessing components/eigen vectors
var = pca.explained_variance_ratio_
print("Amount of variance expained by each PC :\n{}\n".format(var))
print("Sum of variance:{}".format(sum(var)))

cum_var = np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
print("Total variance expained (cumsum) :\n{}".format(cum_var))

# Summary of PCA and Scree/elbow plot (take variables or components from elbow, precipitous to marginal increase)
# Shows proportion of variance explained as a function of principal components(eigen vectors)/predictors:
pca = PCA().fit(x)
plt.plot(cum_var)
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance');
plt.show()
plt.close()

# Including 3rd dataframe for Color reference:
# final_df = pd.concat([principal_comp_df, df_merged[[target]], axis = 1)
output_dir = os.path.expanduser("~/Dropbox/for_genemodels")
output_file = join(output_dir, "PCA_analysis_tf_binding_fraction_reannnotated.txt")
final_df = pd.concat([ principal_comp_df, df_merged[[target]], df_merged[["label"]], df_merged[["tf_name"]], df_merged[["final_ano"]]], axis = 1)

# Rename for later convenience:
final_df["tf_name"] = final_df["tf_name"].str.replace(r"\[FLAG\]", "").str.replace("_human", "")
final_df.rename(columns={0:"Principal_component_1", 1:"Principal_component_2", 2: "Principal_component_3"}, inplace=True)
final_df.to_csv(output_file, sep="\t", header=True, index=False)


###################################
# PCA plot with 3d scatterplot(R) # 
###################################

library(scatterplot3d) # load
library(data.table)
library(dplyr)
library(ggplot2)

output_dir <- "/Users/suryachhetri/Dropbox/for_genemodels/motifs_compinfo_analysis/novel_motif_analysis"
df = fread("/Users/suryachhetri/Dropbox/for_genemodels/PCA_analysis_tf_binding_fraction_reannnotated.txt", sep="\t", header=TRUE)
row.names(df) = df$tf_name
output_file_name = file.path(output_dir,"PCA_clustering_of_TFs.pdf")

# For PCA plot with text:
pdf(output_file_name)

# color = df$label shows color are prelisted in dataframe
s3d = scatterplot3d(df[,1:3], main="Clustering of TFs based on Chromatin Preference", pch = 16, color=df$label) 
s3d.coords = s3d$xyz.convert(df$Principal_component_1, df$Principal_component_2, df$Principal_component_3)
text(s3d.coords$x, s3d.coords$y, labels=row.names(df), pos=4, cex=.18)                  # shrink 
legend("topleft", legend=c("Promoter_like", "Enhancer_like", "Promoter_and_Enhancer_like", "Heterochromatin_and_Repressed", "CTCF_Bound", "Active_Gene_Body"),
col=c("red","orange","#619CFF","gray67","purple3","green"), pch = 16, cex=0.4) 
dev.off()

# For PCA plot without text:
output_file_name = file.path("PCA_analysis_tf_binding_fraction_3D_fig_without_text.pdf")
pdf(output_file_name)
# color = df$label shows color are prelisted in dataframe
scatterplot3d(df[,1:3], main="3-D Spread of TFs based on Chromatin Preference", pch = 16, color=df$label)
legend("right", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
col=c("red","orange","blueviolet","burlywood","green"), pch = 16, cex=1) 
dev.off()


# Plot principal components:
df$final_ano <- factor(df$final_ano, levels=c("Promoter_like", "Enhancer_like", "Promoter_and_Enhancer_like", "Heterochromatin_and_Repressed", "CTCF_Bound", "Active_Gene_Body"))
df$final_ano %>% factor %>% levels
ggplot(df, aes(Principal_component_1, Principal_component_2, color=final_ano, label=tf_name)) + 
geom_point() + geom_text(check_overlap = TRUE, size=1.5, color="black") +
labs(x="Principal Component 1", y="Principal Component 2", color="Annotation") +
ggtitle("2D-PCA")+ #xlim(c(0,1))
scale_color_manual(values=c("red","orange","#619CFF","gray67","purple3","green"))+ 
theme_minimal()

ggsave(file.path(output_dir,"2D_PCA_clustering_of_TFs.pdf"), width = 9, height = 7, dpi=500, device='pdf')
