##########################################
##########################################
# Random Forest Python:
##########################################
##########################################

""" K-fold and Supervised Classification Algorithm 
    for learning models, or generate Predictive models """

import time
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import KFold
import pickle 
from os.path import join
from collections import defaultdict


""" Pre-Processing of dataframe for classifier train """

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")

final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))

# Select "ideas state" of interest:
final_pca_df = final_pca_df.loc[final_pca_df["ideas_state"].isin(["Prom assoc", "Strong Enh"])]

# Initialize Random Forest Classifier:
rf_clf = RandomForestClassifier(n_estimators=500, n_jobs=4, random_state=1)

# 5X KFold and StratifiedKFold Cross-validation
# 10 Repeated 5X KFold Cross-validation (5X StratifiedKFold Cross-validation)
cv = KFold(n_splits=5)
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=1) 
skf = StratifiedKFold(n_splits=5)

# Randomly sample 5000 sets from each state:

predict_dict = defaultdict(list)
acc_score_dict = defaultdict(list)
confusion_mat_dict = defaultdict(list)
acc_score = [] 

start_time = time.time()
for i in range(10):
    idx_iter="Iter {}".format(i)
    print(idx_iter)
    final_pca_df_sample = final_pca_df.groupby('ideas_state').apply(lambda x: x.sample(2500, random_state=i)).reset_index(drop=True)
    final_pca_df_sample["ideas_state"].value_counts()
    #final_pca_df_sample = final_pca_df.copy()

    # Separate features and target/labels from dataframe:
    feature_df = final_pca_df_sample.iloc[:, 6:]
    label_df = final_pca_df_sample["ideas_state"]

    # Tranform features to arrays of features:
    X_feat = feature_df.values 

    # Tranform labels to numpy arrays:
    re_label = label_df.replace({"Prom assoc" : 0, "Strong Enh" : 1})
    y_lab = re_label.values

    i=0
    for train_index, test_index in skf.split(X_feat, y_lab):
        print("\nFold {}".format(i))
        rf_model = rf_clf.fit(X_feat[train_index], y_lab[train_index])
        ypred = rf_model.predict(X_feat[test_index])
        mat = confusion_matrix(y_lab[test_index], ypred)
        acc = accuracy_score(y_lab[test_index], ypred)
        print(acc)
        acc_score.append(acc)
        # print(metrics.classification_report(ypred, y_lab[test_index]))
        predict_dict[idx_iter].append((y_lab[test_index], ypred))
        acc_score_dict[idx_iter].append(acc)
        confusion_mat_dict[idx_iter].append(mat)
        i+=1

end_time = time.time()
print("{} seconds".format(end_time - start_time))

#sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False)


###################################################
##################################################
# Random Forest Train with random feature sampling
##################################################
###################################################

""" K-fold and Supervised Classification Algorithm 
    for learning models, or generate Predictive models """

import time
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import KFold
import pickle 
from os.path import join
from collections import defaultdict
from numpy import random


""" Pre-Processing of dataframe for classifier train """

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")

final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))

# Select "ideas state" of interest:
final_pca_df = final_pca_df.loc[final_pca_df["ideas_state"].isin(["Prom assoc", "Strong Enh"])]

# Initialize Random Forest Classifier:
rf_clf = RandomForestClassifier(n_estimators=500, n_jobs=4, random_state=1)

# 5X KFold and StratifiedKFold Cross-validation
# 10 Repeated 5X KFold Cross-validation (5X StratifiedKFold Cross-validation)
cv = KFold(n_splits=5)
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=1) 
skf = StratifiedKFold(n_splits=5)


predict_dict = defaultdict(list)
acc_score_dict = defaultdict(list)
confusion_mat_dict = defaultdict(list)
#acc_score = [] 

start_time = time.time()
for i in range(3):
    idx_iter="Iter {}".format(i)
    print(idx_iter)
    # Randomly sample 5000(2500 each prom/enhancer) sets from each state:
    final_pca_df_sample = final_pca_df.groupby('ideas_state').apply(lambda x: x.sample(2500, random_state=i)).reset_index(drop=True)

    # Separate features and target/labels from dataframe:
    feature_df = final_pca_df_sample.iloc[:, 6:]
    label_df = final_pca_df_sample["ideas_state"]

    for feat_num in range(1,15,5):
        print("Random Feature Sampling : {}".format(feat_num))
        random.seed(feat_num)
        feature_df_new = feature_df.iloc[:,np.random.randint(0,feature_df.shape[1],size=feat_num)]

        # Tranform features to arrays of features:
        X_feat = feature_df_new.values 

        # Tranform labels to numpy arrays:
        re_label = label_df.replace({"Prom assoc" : 0, "Strong Enh" : 1})
        y_lab = re_label.values

        fold_acc = []
        i=0
        for train_index, test_index in skf.split(X_feat, y_lab):
            print("\nFold {}".format(i))
            rf_model = rf_clf.fit(X_feat[train_index], y_lab[train_index])
            ypred = rf_model.predict(X_feat[test_index])
            mat = confusion_matrix(y_lab[test_index], ypred)
            acc = accuracy_score(y_lab[test_index], ypred)
            print(acc)
            fold_acc.append(acc)
            # print(metrics.classification_report(ypred, y_lab[test_index]))
            predict_dict[idx_iter].append((y_lab[test_index], ypred))
            confusion_mat_dict[idx_iter].append(mat)
            i+=1

        fold_acc_mean = np.mean(fold_acc)
        acc_score_dict[idx_iter].append(fold_acc_mean)

end_time = time.time()
print("{} seconds".format(end_time - start_time))

acc_df = pd.DataFrame.from_dict(acc_score_dict, orient="index")


################################################################
################################################################
# Final Random Forest Train with random feature sampling
# without random iterations; using all datasets for train-test
################################################################
################################################################

""" K-fold and Supervised Classification Algorithm 
    for learning models, or generate Predictive models """

import time
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import KFold
import pickle 
from os.path import join
from collections import defaultdict
from numpy import random
from StringIO import StringIO # 3.0 from io import StringIO

def report_to_df(report):
    report = re.sub(r" +", " ", report).replace("avg / total", "avg/total").replace("\n ", "\n")
    report_df = pd.read_csv(StringIO("Classes" + report), sep=' ', index_col=0)        
    return(report_df)

def metrics_report_to_df(ytrue, ypred):
    precision, recall, fscore, support = metrics.precision_recall_fscore_support(ytrue, ypred)
    classification_report = pd.concat(map(pd.DataFrame, [precision, recall, fscore, support]), axis=1)
    classification_report.columns = ["precision", "recall", "f1-score", "support"] # Add row w "avg/total"
    classification_report.loc['avg/Total', :] = metrics.precision_recall_fscore_support(ytrue, ypred, average='weighted')
    classification_report.loc['avg/Total', 'support'] = classification_report['support'].sum() 
    return(classification_report)

""" Pre-Processing of dataframe for classifier train """

output_dir = expanduser("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm")

final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_2.pkl"))
final_pca_df = pd.read_pickle(join(output_dir, "final_pca_df_10.pkl"))

# Select "ideas state" of interest:
final_pca_df = final_pca_df.loc[final_pca_df["ideas_state"].isin(["Prom assoc", "Strong Enh"])]

# Initialize Random Forest Classifier(# random_state is given for reproducibility)
rf_clf = RandomForestClassifier(n_estimators=500, n_jobs=4, random_state=1) 

# 5X KFold and StratifiedKFold Cross-validation
# 10 Repeated 5X KFold Cross-validation (5X StratifiedKFold Cross-validation)
cv = KFold(n_splits=5)
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=1) 
skf = StratifiedKFold(n_splits=5)

ypred_dict = defaultdict(list)
yscore_dict = defaultdict(list)
acc_score_dict = defaultdict(list)
confusion_mat_dict = defaultdict(list)
clf_report_dict = defaultdict(list)

# Separate features and target/labels from dataframe:
feature_df = final_pca_df.iloc[:, 6:]
label_df = final_pca_df["ideas_state"]

start_time = time.time()
for feat_num in range(5,205,5):
    print("Random Feature Sampling : {}".format(feat_num))
    random.seed(feat_num+1) # random_state seed
    feature_df_new = feature_df.iloc[:,np.random.randint(0,feature_df.shape[1],size=feat_num)]

    # Tranform features to arrays of features:
    X_feat = feature_df_new.values 

    # Tranform labels to numpy arrays:
    re_label = label_df.replace({"Prom assoc" : 0, "Strong Enh" : 1})
    y_lab = re_label.values

    i=0
    for train_index, test_index in rskf.split(X_feat, y_lab):
        print("\nFold {}".format(i))
        rf_model = rf_clf.fit(X_feat[train_index], y_lab[train_index])
        yscore = rf_model.predict_proba(X_feat[test_index])[:,1] # Considering 1(Enh) as positive label
        ypred = rf_model.predict(X_feat[test_index])
        clf_report = metrics.classification_report(y_lab[test_index], ypred) # Esp. useful for multiclass clf
        clf_report_df = report_to_df(clf_report).reset_index()
        mat = confusion_matrix(y_lab[test_index], ypred)
        acc = accuracy_score(y_lab[test_index], ypred)
        print(acc)
        acc_score_dict[feat_num].append(acc)
        clf_report_dict[feat_num].append(clf_report_df)
        confusion_mat_dict[feat_num].append(mat)
        ypred_dict[feat_num].append((y_lab[test_index], ypred))
        yscore_dict[feat_num].append((y_lab[test_index], yscore))
        #clf_report_df = pd.read_csv(StringIO(clf_report), sep="\t")
        i+=1

end_time = time.time()
print("{} seconds".format(end_time - start_time))

acc_df = pd.DataFrame.from_dict(acc_score_dict, orient="index").sort_index()
acc_df.columns = ["fold " + str(each) for each in range(10)] # for 5X crossfold-val
acc_df = acc_df.reset_index().rename(columns={"index":"feat_index"})

# Accuracy mean across rows :(sorted by features)
acc_mean = acc_df.iloc[:,1:].apply(lambda row: row.mean(), axis=1)
final_acc_df = pd.concat([acc_df, acc_mean], axis=1).rename(columns={0: "acc_mean"})
print(acc_mean); print(final_acc_df)
  
# Mock dataset for plot to find min of row and max of row (across fold-change for arrow-point plot)
row_min = final_acc_df.iloc[:,1:-1].apply(lambda row: row.min()-0.4,axis=1)
row_max = final_acc_df.iloc[:,1:-1].apply(lambda row: row.max()+0.1,axis=1)


####################################

from plotnine import *

p = (ggplot(final_acc_df, aes(y="acc_mean", x="feat_index")) +
        geom_point(color="red") + 
        geom_point(aes(y=acc_mean-0.2), color="blue")+ 
        geom_segment(aes(y=row_min, x="feat_index", yend=row_max, xend="feat_index"), color="grey") +
        ylim(0,1))
# Find out the order of legends using manual factoring:
pd.Categorical(["geom_2", "geom_1", "geom_3"])

# For color change and label change manual 
scale_colour_manual(
    name="Annotation", 
    labels = ["mean", "mean_val-0.2"], 
    values = ["black", "blue"])

#####################################


# Confusion matrix and classfication report useful for multiclass-multilabel problems:
pr_recall_f1_list = []
for key, value in clf_report_dict.items():
    report_df = pd.concat(clf_report_dict[key], ignore_index=True)
    report_df["feat_index"] = key
    pr_recall_f1_list.append(report_df) 

metrics_report_df = pd.concat(pr_recall_f1_list, ignore_index=True)
metrics_report_sort = metrics_report_df.sort_values(["feat_index", "Classes"])
pr_report = metrics_report_sort.groupby(["feat_index", "Classes"]).apply(lambda x: x["precision"].mean()).reset_index(name="precision")
recall_report = metrics_report_sort.groupby(["feat_index", "Classes"]).apply(lambda x: x["recall"].mean()).reset_index(name="recall")
f1_score_report = metrics_report_sort.groupby(["feat_index", "Classes"]).apply(lambda x: x["f1-score"].mean()).reset_index(name="f1-score")
support_report = metrics_report_sort.groupby(["feat_index", "Classes"]).apply(lambda x: x["support"].sum()).reset_index(name="support")
pr_recall_f1_report = reduce(lambda left,right: pd.merge(left,right), [pr_report, recall_report,f1_score_report])


#### ROC, PR, AUROC, AUPR for each fold for any given feature: #####

auroc_dict = defaultdict(list)
aupr_dict = defaultdict(list)
roc_curve_dict = defaultdict(list)
pr_curve_dict = defaultdict(list)

# Calculate mean Area-under-curve for ROC/PR and store as df:
for key, value in yscore_dict.items():
    import copy # make deep copy while coypying list, to avoid updating to dict itself - if any
    yscore_pred_list = copy.deepcopy(yscore_dict[key])
    # Find y_true(true label) and y_score (predicted score)
    for each_array in yscore_pred_list:
        y_true, y_scores = each_array
        # print y_true, y_scores

        """ For Area Curve using direct function from sklearn: """
        # pr_areacurve = metrics.average_precision_score(y_true, y_scores)
        # roc_areacurve = metrics.roc_auc_score(y_true, y_scores)
        # aupr_dict[key].append(pr_areacurve)
        # auroc_dict[key].append(roc_areacurve)

        """ For plotting Curve and AUC by calculating precision, recall, tpr, fpr (from scratch):"""
        precision, recall, thresholds = metrics.precision_recall_curve(y_true, y_scores) # y_scores[:,1] for gkm-SVM
        pr_curve_df = pd.DataFrame({"precision" : precision, "recall" : recall}) 
        pr_areacurve = auc(recall, precision)
        aupr_dict[key].append(pr_areacurve)
        pr_curve_dict[key].append(pr_curve_df)

        fpr, tpr, thresholds = metrics.roc_curve(y_true, y_scores) # y_scores[:,1] for gkm-SVM
        roc_curve_df = pd.DataFrame({"true_pos_rate" : tpr, "false_pos_rate" : fpr}) 
        roc_areacurve = auc(fpr, tpr)
        auroc_dict[key].append(roc_areacurve)
        roc_curve_dict[key].append(roc_curve_df)

auroc_df = pd.DataFrame.from_dict(auroc_dict, orient="index").sort_index()
auroc_df.columns = ["fold " + str(each) for each in range(10)] # for 5X crossfold-val
auroc_df = auroc_df.reset_index().rename(columns={"index":"feat_index"})
auroc_mean = auroc_df.iloc[:,1:].apply(lambda row: row.mean(), axis=1)
auroc_mean_df = pd.DataFrame({"feat_index" : auroc_df.iloc[:,0], "auroc_mean": auroc_mean})

aupr_df = pd.DataFrame.from_dict(aupr_dict, orient="index").sort_index()
aupr_df.columns = ["fold " + str(each) for each in range(10)] # for 5X crossfold-val
aupr_df = aupr_df.reset_index().rename(columns={"index":"feat_index"})
aupr_mean = aupr_df.iloc[:,1:].apply(lambda row: row.mean(), axis=1)
aupr_mean_df = pd.DataFrame({"feat_index" : aupr_df.iloc[:,0], "auroc_mean": aupr_mean})

# Calculate Curve for ROC and PR for latter plotting - stored as df:
roc_curve_df_list = []
for key, value in roc_curve_dict.items():
    import copy # make deep copy while coypying list, to avoid updating to dict itself
    df_list = copy.deepcopy(roc_curve_dict[key]) 
    for idx, each_roc_df in enumerate(df_list):
        each_roc_df["feat_index"] = key
        each_roc_df["fold"] = "fold " + str(idx)
        roc_curve_df_list.append(each_roc_df)
        #print (idx, each)

# ROC Curve dataframe plotting data points for each feature and folds:
roc_curve_combined_df = pd.concat(roc_curve_df_list, ignore_index=True)
roc_curve_sorted_df = roc_curve_combined_df.sort_values(["feat_index", "fold"]).reset_index(drop=True)
roc_curve_tfcount_200 = roc_curve_sorted_df.loc[roc_curve_sorted_df["feat_index"] == 200]

pr_curve_df_list = []
for key, value in pr_curve_dict.items():
    import copy # make deep copy while coypying list, to avoid updating to dict itself
    df_list = copy.deepcopy(pr_curve_dict[key]) 
    for idx, each_pr_df in enumerate(df_list):
        each_pr_df["feat_index"] = key
        each_pr_df["fold"] = "fold " + str(idx)
        pr_curve_df_list.append(each_pr_df)
        #print (idx, each)

# PR Curve dataframe plotting data points for each feature and folds:
pr_curve_combined_df = pd.concat(pr_curve_df_list, ignore_index=True)
pr_curve_sorted_df = pr_curve_combined_df.sort_values(["feat_index", "fold"]).reset_index(drop=True)
pr_curve_tfcount_200 = pr_curve_sorted_df.loc[pr_curve_sorted_df["feat_index"] == 200]


#################################################################
#################################################################


""" Direct cross-validation library approach"""

# from sklearn import datasets
# from sklearn.cross_validation import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score, cross_val_predict
# from sklearn.cross_validation import cross_val_score, cross_val_predict
# from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn import metrics
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

# Initialize Random Forest Classifier:
rf_model = RandomForestClassifier(n_estimators=500)

# Define functions:
def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]
def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]

# Scoring metrics for model evaluation: 
scoring = {'accuracy': make_scorer(accuracy_score),
            'prec': 'precision', 'rec' : 'recall', 
            'tp' : make_scorer(tp), 'tn' : make_scorer(tn),
            'fp' : make_scorer(fp), 'fn' : make_scorer(fn),
            'AUC': 'roc_auc', "PRAUC": 'average_precision'}

# Using Cross-validate function (for easy scoring dict arrays)
cv_results = cross_validate(rf_model.fit(X_feat, y_lab), X_feat, y_lab, cv=3, 
            scoring=scoring,return_train_score=False) 

#################

# Using cross_val_predict to access all data point performance
y_true = y_lab.copy()
y_pred = cross_val_predict(rf_model, X_feat, y_lab, cv=5)
accuracy = accuracy_score(y_pred.astype(int), y_lab.astype(int)); print(accuracy)
mat = confusion_matrix(y_true, y_pred)

# Using cross_val_score to access all data point performance
from sklearn.cross_validation import KFold
kf = KFold(len(y_lab), n_splits=5, random_state=0)
cval_score = cross_val_score(rf_model, X_feat, y_lab, cv=kf) #customised
cval_score = cross_val_score(rf_model, X_feat, y_lab, cv=5) #default stratifiedKfold

from sklearn.cross_validation import StratifiedKFold
skf = StratifiedKFold(y_lab, n_folds=5, random_state=0)
cval_score = cross_val_score(rf_model, X_feat, y_lab, cv=skf) #customised
cval_score = cross_val_score(rf_model, X_feat, y_lab, cv=5) #default stratifiedKfold

# Using cross_val_predict (with proba) for ROC and PRAUC to access all data point performance
y_true = y_lab.copy()
y_scores = cross_val_predict(rf_model, X_feat, y_lab, method='predict_proba', cv=skf)

# AUC and PRAUC
from sklearn.metrics import precision_recall_curve, roc_curve, auc
precision, recall, thresholds = precision_recall_curve(y_true, y_scores[:,1])
pr_areacurve = auc(recall, precision)

fpr, tpr, thresholds = roc_curve(y_true, y_scores[:,1])
roc_areacurve = auc(fpr, tpr)

print('Area Under PR Curve(AP): {0:0.4f}'.format(pr_areacurve))
print('Area Under ROC Curve(AP): {0:0.4f}'.format(roc_areacurve))

#pr_areacurve = average_precision_score(y_true, y_scores) 
#roc_areacurve = roc_auc_score(y_true, y_scores)


###############################################################
###############################################################
### Misc:

idx = np.arange(0, len(y))
for j in np.random.randint(0, high=10000, size=10):
    np.random.shuffle(idx)
    kf = KFold(n=len(y), n_folds=10, random_state=j)

    for i, (train, test) in enumerate(kf):
        model = LogisticRegression().fit(X[idx][train], y[idx][train])
        y_score = model.predict_proba(X[idx][test])
        fpr, tpr, _ = roc_curve(y[idx][test], y_score[:, 1])
        plt.plot(fpr, tpr, 'b', alpha=0.05)
        tpr = interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)


# Data extraction approach for feature of interest with AUROC, AUPR:
y_true, y_pred = ypred_dict[5][0] # here 5 represents key (with 5 features) [fold 1]
y_true, y_score = yscore_dict[20][1] # here 20 represents key (with 20 features) [fold 2]

# AUC and PRAUC (if needed for each folds-crossval)
from sklearn.metrics import precision_recall_curve, roc_curve, auc

# For class-imbalance or data-imbalance case:
precision, recall, thresholds = precision_recall_curve(y_true, y_scores[:,1])
pr_areacurve = auc(recall, precision)
print('Area Under PR Curve(AP): {0:0.4f}'.format(pr_areacurve))

# For class-balance or data-balance case:
fpr, tpr, thresholds = roc_curve(y_true, y_scores[:,1])
roc_areacurve = auc(fpr, tpr)
print('Area Under ROC Curve(AP): {0:0.4f}'.format(roc_areacurve))

##### Lollipop plot (or arrow point plot) with dot as primary reference and 
##### removing thick bar with thin bar to replace bulky barplots

######## Using R ###############

library(dplyr)          # for data manipulation
library(tidyr)          # for data tidying
library(ggplot2)        # for generating the visualizations

# Test dataset:
head(midwest)

# Mock data from R:
ohio_top25 <- midwest %>%
        filter(state == "OH") %>%
        select(county, percollege) %>%
        arrange(desc(percollege)) %>%
        top_n(25) %>%
        arrange(percollege) %>%
        mutate(county = factor(county, levels = .$county))

# Dot plot
ggplot(ohio_top25, aes(percollege, county)) +
        geom_point()

# Simple Lollipop plot (single point):
ggplot(ohio_top25, aes(percollege, county)) +
        geom_segment(aes(x = 0, y = county, xend = ohio_top25$percollege, yend = county), color = "grey50") +
        xlim(0,50) + geom_point(col="blue") # (puts point on X with respect to y(i.e county here))

# cols <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")
# geom_errorbar(aes(ymin=Avg-AvgSE, ymax=Avg+AvgSE))

p <- ggplot(ohio_top25, aes(percollege, county)) +
        geom_segment(aes(x = min(ohio_top25$percollege), y = county, xend = max(ohio_top25$percollege), yend = county), color = "grey50") +
        xlim(0,50) + 
        # Manually adding color to ggplot:
        geom_point(aes(ohio_top25$percollege+5, county, color="(meanval+5)"), pch=16) +
        geom_point(aes(median(ohio_top25$percollege), county, color="(median)"), pch=8) + 
        theme(axis.text.x=element_text(angle=30, hjust=1))
        #scale_colour_manual(name="Error Bars",values=cols) + scale_fill_manual(name="Bar",values=cols)

p + scale_colour_manual(name="Annotation", values = c("black", "red")) + 
  guides(colour = guide_legend(override.aes = list(shape = c(16, 8))))

# Implementing plotting with real dataset (Borrowed from pandas df) for R ggplot:
final_acc_df <- fread(file.path(output_dir, "RF_accuracy_mean.txt"), sep="\t") %>% as.data.frame
df_filt <- final_acc_df[,3:ncol(df)-1]
row_min <- apply(df_filt,1,min)-0.2
row_max <- apply(df_filt,1,max)+0.1

p <- ggplot(final_acc_df, aes(y=acc_mean, x=feat_index)) +
        geom_point(aes(color="geom_1"),pch=8)  + 
        geom_segment(aes(y=row_min, x=feat_index, yend=row_max, xend=feat_index), color="black") + 
        geom_point(aes(y=acc_mean+0.1, x=feat_index, color="geom_2"), pch=16) +
        geom_point(aes(y=acc_mean-0.1, x=feat_index, color="geom_3"), pch=16) +
        ylim(0,1)

# Find out the order of legends using manual factoring:
factor(c("geom_2", "geom_1", "geom_3"))
# output = Levels: geom_1, geom_2, geom_3

p + scale_colour_manual(
    name="Annotation", 
    labels = c("mean", "mean_val+0.1", "mean_val-0.1"), 
    values = c("black", "blue", "darkred")) + 
    guides(colour = guide_legend(override.aes = list(shape = c(8, 16, 16))))
    # + theme(legend.position = "top")





