library(ranger)

#Make binary matrix of loci by factor bound
HighBound<-read.table("~/Dropbox/for_ryne/random_forest/MergedBed_2TFcutoff_IDEAS.txt",sep="\t",header=T,stringsAsFactors = F)
Files<-list.files("~/Dropbox/kmers_svm/UniqueENCODETFs_renamed/")
#Files<-list.files("~/Surya_HEPG2_Beds/")

HotFrame<-matrix(nrow=nrow(HighBound),ncol=length(Files))
colnames(HotFrame)<-gsub(".bed","",Files)
rownames(HotFrame)<-paste(HighBound[,1],HighBound[,2],HighBound[,3],sep="_")
rownames(HighBound)<-paste(HighBound[,1],HighBound[,2],HighBound[,3],sep="_")

# Forming binary matrix of 0/1 using Hotsites:
for(i in colnames(HotFrame)){
  print(i)
  HotFrame[grep(i,HighBound[,5]),i]<-1
}
HotFrame[which(is.na(HotFrame))]<-0


#Filter Random forest input to desirable HOTness
HighBound_Filt<-HighBound[which(HighBound[,6]>10),]
HotFrame_Filt<-HotFrame[which(HighBound[,6]>10),]

#Run Random forest loop
Master_Perf<-list()
for(i in 1:100){
  print(i)
  set.seed(i)
  PromSamp<-sample(grep("prom_assoc",HighBound_Filt[,7]),1000)
  set.seed(i)
  EnhSamp<-sample(grep("sEnh_assoc",HighBound_Filt[,7]),1000)
  Perf<-c()
  for(DropNum in seq(0,208,5)){
    #print(DropNum)
    toTrain<-HotFrame_Filt[c(PromSamp,EnhSamp),]
    set.seed(i)
    toTrain<-toTrain[,sample(1:ncol(toTrain),ncol(toTrain)-DropNum)]
    #rMod<-randomForest(x=toTrain,y=as.factor(c(rep(0,1000),rep(1,1000))))
    
    y=as.factor(c(rep(0,1000),rep(1,1000)))
    rMod<-ranger(y~., data=data.frame(toTrain), importance="impurity")
    Perf[paste0(DropNum)]<-mean(data.frame(rMod$confusion)[2:3,3])
  }
  Master_Perf[[i]] <- Perf
}

Master_Perf_Frame<-do.call(rbind, Master_Perf)
plot(apply(Master_Perf_Frame,2,median)/10, xaxt="n", ylab="Percent Missclassified",ylim=c(0,50),xlab="Number of Factors Excluded",pch=16, col="red")
arrows(x0 = 1:42,x1 = 1:42,y0 = apply(Master_Perf_Frame,2,min)/10,y1 = apply(Master_Perf_Frame,2,max)/10,angle = 90,code = 3,length = 0.05)
axis(1, labels=seq(0,207,5)[seq(1,42,5)],at=seq(1,42,5))


##############################################################
##############################################################
##############################################################
# Percent misclassified approach, using train (OOB error)

library(ranger)
library(data.table)
library(dplyr)
library(stringr)

# Make binary matrix of loci by factor bound
HighBound_Filt<-read.table("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hotsite_10_TF_random_forest_train.txt",sep="\t",header=T, stringsAsFactors = F)
HighBound_Filt<-read.table("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Binary_Hotsite_10_TF_random_forest_train.txt",sep="\t",header=T, stringsAsFactors = F)
HotFrame_Filt <- HighBound_Filt[,7:ncol(HighBound_Filt)]

# Rename columns for DBF filtering:
names(HotFrame_Filt) <- names(HotFrame_Filt) %>% str_replace("TF_", "") %>% str_replace("\\.", "[") %>% str_replace("\\.", "]")

# Load TF annoation file:
anno_df <- fread("~/Dropbox/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))
dbf_anno_df <- anno_df %>% filter(Category=="DBF")

# Include only DBF's or variables v1, v2, v3:
myvars <- names(HotFrame_Filt) %in% dbf_anno_df$Target
dbf_HotFrame_Filt <- HotFrame_Filt[myvars]

# Run Random forest loop
Master_Perf_motif<-list()
for(i in 1:100){
  print(i)
  set.seed(i)
  Prom_Samp<-sample(grep("Prom assoc",HighBound_Filt[,6]),1000)
  set.seed(i)
  Enh_Samp<-sample(grep("Strong Enh",HighBound_Filt[,6]),1000)
  
  Perf<-c()
  exclude_factors_int <- seq(0,165,5)
  for(DropNum in exclude_factors_int){
    #print(DropNum)
    toTrain<-dbf_HotFrame_Filt[c(Prom_Samp, Enh_Samp),]
    set.seed(i)
    toTrain<-toTrain[,sample(1:ncol(toTrain),ncol(toTrain)-DropNum)]
    #rMod<-randomForest(x=toTrain,y=as.factor(c(rep(0,1000),rep(1,1000))))
    
    y=as.factor(c(rep(0,1000),rep(1,1000)))
    rMod<-ranger(y~., data=data.frame(toTrain), importance="impurity")
    Perf[paste0(DropNum)]<-mean(data.frame(rMod$confusion)[2:3,3])
  }
  Master_Perf_motif[[i]]<-Perf
}

# Master_Perf_motif_Frame<-do.call(rbind, Master_Perf_motif)
pdf("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Random_forest_trained_motif_plot_1_binary.pdf")
plot(apply(Master_Perf_motif_Frame,2,median)/10, xaxt="n", ylab="Element Missclassified (%)",ylim=c(0,50),xlab="Number of Factors Excluded",pch=16, col="red", main="Random Forest Trained on Factors Binding Matrix")
arrows(x0 = 1:length(exclude_factors_int),x1 = 1:length(exclude_factors_int),y0 = apply(Master_Perf_motif_Frame,2,min)/10,y1 = apply(Master_Perf_motif_Frame,2,max)/10,angle = 90,code = 3,length = 0.05)
axis(1, labels=exclude_factors_int[seq(1,length(exclude_factors_int),5)],at=seq(1,length(exclude_factors_int),5))
abline(h=15, col="blue", lwd=1, lty=2)
abline(v=11, col="blue", lwd=1, lty=2)
dev.off()

##############################################################
# Percent accurracy approach, using train (OOB error)

library(ranger)
library(data.table)
library(dplyr)
library(stringr)

# Make binary matrix of loci by factor bound
HighBound_Filt<-read.table("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Hotsite_10_TF_random_forest_train.txt",sep="\t",header=T, stringsAsFactors = F)
#HighBound_Filt<-read.table("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Binary_Hotsite_10_TF_random_forest_train.txt",sep="\t",header=T, stringsAsFactors = F)
HotFrame_Filt <- HighBound_Filt[,7:ncol(HighBound_Filt)]

# Rename columns for DBF filtering:
names(HotFrame_Filt) <- names(HotFrame_Filt) %>% str_replace("TF_", "") %>% str_replace("\\.", "[") %>% str_replace("\\.", "]")

# Load TF annoation file:
anno_df <- fread("~/Dropbox/misc/TFs_Annotation_file.txt", sep="\t")
anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))
dbf_anno_df <- anno_df %>% filter(Category=="DBF")

# Include only DBF's or variables v1, v2, v3:
myvars <- names(HotFrame_Filt) %in% dbf_anno_df$Target
dbf_HotFrame_Filt <- HotFrame_Filt[myvars]

# Run Random forest loop
Master_Perf_motif<-list()
for(i in 1:100){
  print(i)
  set.seed(i)
  Prom_Samp<-sample(grep("Prom assoc",HighBound_Filt[,6]),1000)
  set.seed(i)
  Enh_Samp<-sample(grep("Strong Enh",HighBound_Filt[,6]),1000)
  
  Perf<-c()
  include_factors_int <- seq(1,170,5)
  for(includeNum in include_factors_int){
    #print(includeNum)
    toTrain<-dbf_HotFrame_Filt[c(Prom_Samp, Enh_Samp),]
    set.seed(i)
    toTrain<-toTrain[,sample(1:ncol(toTrain), includeNum)]
    #rMod<-randomForest(x=toTrain,y=as.factor(c(rep(0,1000),rep(1,1000))))
    
    y=as.factor(c(rep(0,1000),rep(1,1000)))
    rMod<-ranger(y~., data=data.frame(toTrain), importance="impurity")
    Perf[paste0(includeNum)]<-mean(data.frame(rMod$confusion)[c(1,4),3])
  }
  Master_Perf_motif[[i]]<-Perf
}

Master_Perf_motif_Frame_include<-do.call(rbind, Master_Perf_motif)
pdf("~/Dropbox/for_genemodels/motifs_compinfo_analysis/kmer_svm/Random_forest_trained_motif_plot_2_binary.pdf")
plot(apply(Master_Perf_motif_Frame_include,2,median)/10, xaxt="n", ylab="Element classification accurracy (%)",ylim=c(0,100),xlab="Number of Factors Included",pch=16, col="red", main="Random Forest Trained on Factors Binding Matrix")
arrows(x0 = 1:length(include_factors_int),x1 = 1:length(include_factors_int),y0 = apply(Master_Perf_motif_Frame_include,2,min)/10,y1 = apply(Master_Perf_motif_Frame_include,2,max)/10,angle = 90,code = 3,length = 0.05)
axis(1, labels=c(1,seq(5,170,5))[seq(1,length(include_factors_int),5)],at=seq(1,length(include_factors_int),5))
abline(h=80, col="blue", lwd=1, lty=2)
abline(v=9,  col="blue", lwd=1, lty=2)
dev.off()

