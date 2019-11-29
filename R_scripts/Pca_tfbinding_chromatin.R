
###################################
# PCA plot with 3d scatterplot(R) # 
###################################


library(scatterplot3d) # load
library(data.table)

df = fread("PCA_analysis_tf_binding_fraction.txt", sep="\t", header=TRUE)
row.names(df) = df$tf_name
output_file_name = file.path("PCA_analysis_tf_binding_fraction_3D_wNames_fig.pdf")

# For PCA plot with text:
pdf(output_file_name)

# color = df$label shows color are prelisted in dataframe
s3d = scatterplot3d(df[,1:3], main="3-D Spread of TFs based on Chromatin Preference", pch = 16, color=df$label) 

s3d.coords = s3d$xyz.convert(df$Principal_component_1, df$Principal_component_2, df$Principal_component_3)
text(s3d.coords$x, s3d.coords$y, labels=row.names(df), pos=4, cex=.18)                  # shrink 

legend("right", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
col=c("red","orange","blueviolet","burlywood","green"), pch = 16, cex=1) 
dev.off()


# For PCA plot without text:
output_file_name = file.path("PCA_analysis_tf_binding_fraction_3D_fig.pdf")
pdf(output_file_name)

# color = df$label shows color are prelisted in dataframe
scatterplot3d(df[,1:3], main="3-D Spread of TFs based on Chromatin Preference", pch = 16, color=df$label)

legend("right", legend=c("Prom assoc","Enh assoc","CTCF assoc","Multiple assoc","Others"),
col=c("red","orange","blueviolet","burlywood","green"), pch = 16, cex=1) 
dev.off()


# Scatter plot to project 36 dimensions to 2 dim
# plt.scatter(final_df.iloc[:, 0].values, final_df.iloc[:, 1].values,
#             c=final_df.loc[:,target].values, edgecolor='none', alpha=0.5,
#             cmap=plt.cm.get_cmap('spectral', 5))

# plt.xlabel('component 1')
# plt.ylabel('component 2')
# plt.colorbar();

# targets = ['Iris-setosa', 'Iris-versicolor', 'Iris-virginica']
# df = finalDf.loc[finalDf["target"].isin(targets)]

