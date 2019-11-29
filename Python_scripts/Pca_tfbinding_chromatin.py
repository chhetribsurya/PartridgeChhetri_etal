
##########################################################
# PCA analysis based on TFs binding at chromatin states
##########################################################


import pandas as pd, numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import matplotlib; import matplotlib.pyplot as plt
import os; from os.path import join, expanduser
from ggplot import *
from plotnine import *

input_file = expanduser("~/Dropbox/for_genemodels/ideas_piestate_dist_for_heatmap.txt")
df_orig = pd.read_csv(input_file, sep="\t") # original value containing file
df_new = pd.read_csv("~/Dropbox/TFs_Annotation_file.txt", sep="\t") # TF annotation file
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
output_file = join(output_dir, "PCA_analysis_tf_binding_fraction.txt")
final_df = pd.concat([ principal_comp_df, df_merged[[target]], df_merged[["label"]], df_merged[["tf_name"]] ], axis = 1)

# Rename for later convenience:
final_df["tf_name"] = final_df["tf_name"].str.replace("\[FLAG\]", "").str.replace("_human", "")
final_df.rename(columns={0:"Principal_component_1", 1:"Principal_component_2", 2: "Principal_component_3"}, inplace=True)
final_df.to_csv(output_file, sep="\t", header=True, index=False)

# Order the factor/categorical variable to color legend accordingly:
final_df["annotation"] = pd.Categorical(final_df["annotation"], categories=["Ideas Prom", "Ideas Enh,", "Ideas CTCF", "Ideas multiple", "Ideas Other"], ordered=True)

# Plot principal components:
plot = (ggplot(final_df) + 
		aes("Principal_component_1", "Principal_component_2", fill=target)+ 
		geom_point() +
		ggtitle("2 Component PCA") +
		# theme_bw() +
		scale_x_continuous(name="Principal Component 1") + scale_y_continuous(name="Principal Component 2") +
		# xlab("Principal Component 1") + ylab("Principal Component 2") +
		scale_fill_manual(name="IDEAS Anno",values=["red","orange","blueviolet","burlywood","green"]) +
		# guides(fill=guide_legend(title="IDEAS Anno"))
		theme(
		axis_title_y = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		axis_title_x = element_text(size=12), #theme(axis_text_x=element_text(angle=45))
		plot_title = element_text(size=14, face="bold"),
		legend_title = element_text(size=8, face="bold")
		# axis_text_y = element_text(size=1.3),
		# axis_text_y = element_text(size=1.3) 
		)
	)

ggsave(plot,"PCA_analysis_tf_binding_fraction_2D.pdf", width=6, height=6)

