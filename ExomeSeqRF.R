rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Esome Sequencing Donor/Recipient pairs
###
### CITATION: 
###
### PROCESS: Random Forest analysis
###           
### DESCRIP: ExomeSeq Data from Donors/Recipients pairs 
###           27 pairs + 1 pair with two donors
###         
###
### Author: Silvia Pineda
### Date: June, 2016
############################################################################################

library("ggplot2")
library("gridExtra")
library("cowplot")
library("randomForest")
library("RColorBrewer")
library("ROCR")
library("party")
library("pROC")

#Work Directory
setwd("/Users/Pinedasans/Catalyst/Results/")
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq_Diff_demo.Rdata")
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq/ExomeSeqVCF_SNPs.Rdata")
demographics<-read.table("/Users/Pinedasans/Catalyst/Data/Demographics.txt",sep="\t",header=T)
non.list<- seq(1,56,2) ##Donors

#########################
##### Random Forest #### 
########################
#exome_variants_diff_imputed<-rfImpute(data.frame(exome_variants_diff2),demographics$phenotype[non.list])
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq/exome_variants_diff_imputed.Rdata")

##Prepare data to run VSURF in the server
exome_variants_diff_complete<-exome_variants_diff2[ , ! apply( exome_variants_diff2, 2 , function(x) any(is.na(x)) ) ] #450,981
save(exome_variants_diff_complete,demographics,file="/Users/Pinedasans/Data/Catalyst/DataRF.Rdata")

####################
### 1a. Apply random forest to find the variants selected by VSURF using three categories (AMR, CMR, NoRej) with the imputed data
load("/Users/Pinedasans/Catalyst/Results/ResultsExomeRFimputed.Rdata")

# Create model with the variants found with the VSURF
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",3,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,1))
dataset<-data.frame(exome_variants_diff_imputed[,fit$varselect.interp])

set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-demographics$phenotype[non.list][-i]
  testDataClass<-demographics$phenotype[non.list][i]
  
  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData)) # Prediction
}
result_roc <- multiclass.roc(endpoint, result_predicted)$auc

exome_variants_rf_selected<-exome_variants_diff_imputed[,fit$varselect.interp]
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected),proximity=TRUE, keep.forest=T,ntree=100)

##plot MDS
tiff("/Users/Pinedasans/Catalyst/Article/MDSplot_discovery.tiff", width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "Discovery",cex=1.6)
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20,pt.cex=1.6 ,cex=1.2)
dev.off()

###Apply logistic regression to find only those that are risk 
exome_variants_rf_selected<-exome_variants_diff_complete[,fit$varselect.interp]
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",1,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,3))
OR<-NULL
p.value.OR<-NULL
for (i in 1:ncol(exome_variants_rf_selected)){
  print(i)
  model<-glm(exome_variants_rf_selected[,i]~endpoint,family = "binomial")
  #m<-polr(factor(endpoint) ~ exome.variants.sign[,i],Hess=TRUE)
  OR[i]<-exp(coef(model))[2]
  p.value.OR[i]<-coefficients(summary(model))[2,4]
}

Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_rf_selected)){
  Diff.AMR[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_rf_selected),df_joint_qc$snp_id)
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,OR,p.value.OR,rf_output$importance)
write.table(df_joint_sign,file="/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",sep="\t",row.names = F)

####################
### 1b. Apply random forest to find the variants selected by VSURF using three categories (AMR, CMR, NoRej)
load("/Users/Pinedasans/Catalyst/Results/ResultsRF.Rdata")

# Create model with the variants found with the VSURF
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",3,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,1))
dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])

set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-demographics$phenotype[non.list][-i]
  testDataClass<-demographics$phenotype[non.list][i]
  
  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData)) # Prediction
}
result_roc <- multiclass.roc(endpoint, result_predicted)$auc

exome_variants_rf_selected<-exome_variants_diff_complete[,fit$varselect.interp]
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected),proximity=TRUE, keep.forest=T,ntree=100)

###To find the trees within random forest
split_var<-list()
for (i in 1:100){
  tree <- getTree(rf_output_total, i, labelVar=TRUE)
  split_var[[i]]<-tree$`split var`
}
tree <- getTree(rf_output_total, 1, labelVar=TRUE)
plot.getTree(rf_output_total, k=8, depth=4)

##plot MDS
tiff("/Users/Pinedasans/Catalyst/Article/MDSplot_discovery.tiff", width = 6, height = 6, units = 'in', res = 300, compression = 'lzw')
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "Discovery",cex=1.6)
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20,pt.cex=1.6 ,cex=1.2)
dev.off()

###Apply logistic regression to find only those that are risk 
exome_variants_rf_selected<-exome_variants_diff_complete[,fit$varselect.interp]
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",1,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,3))
OR<-NULL
p.value.OR<-NULL
for (i in 1:ncol(exome_variants_rf_selected)){
  print(i)
  model<-glm(exome_variants_rf_selected[,i]~endpoint,family = "binomial")
  #m<-polr(factor(endpoint) ~ exome.variants.sign[,i],Hess=TRUE)
  OR[i]<-exp(coef(model))[2]
  p.value.OR[i]<-coefficients(summary(model))[2,4]
}

Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_rf_selected)){
  Diff.AMR[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_rf_selected[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_rf_selected),df_joint_qc$snp_id)
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,OR,p.value.OR,rf_output$importance)
write.table(df_joint_sign,file="/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",sep="\t",row.names = F)


####################################
####Heatmap with the results#######
##################################
library("pheatmap")
race <- ifelse(demographics$mismatch[non.list] == "0","Non race-mismatch", "race-mismatch")
related <- ifelse(demographics$REL2[non.list] == "0","Non related", "related")

annotation_row = data.frame(
  phenotype = demographics$phenotype[non.list],
  race = race,
  related = related)
rownames(annotation_row)<-variant_list$Id2
rownames(exome_variants_rf_selected)<-variant_list$Id2
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
fill=brewer.pal(3,"Set1")
ann_colors = list (phenotype = c(AMR = fill[1], CMR = fill[2], "No-REJ" = fill[3]),
                   related = c("Non related" = "yellow", related = "lightblue"))
pheatmap(exome_variants_rf_selected ,cluster_rows = T,color = colorRampPalette(c("white", "red"))(50),cex=1.1,
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors,border_color=F)


###plot the tree
options(repos='http://cran.rstudio.org')
have.packages <- installed.packages()
cran.packages <- c('devtools','plotrix','randomForest','tree')
to.install <- setdiff(cran.packages, have.packages[,1])
if(length(to.install)>0) install.packages(to.install)

library(devtools)
if(!('reprtree' %in% installed.packages())){
  install_github('araastat/reprtree')
}
for(p in c(cran.packages, 'reprtree')) eval(substitute(library(pkg), list(pkg=p)))

reprtree:::plot.getTree(rf_output_total)
tree <- getTree(rf_output_total, k=1, labelVar=TRUE)
realtree <- reprtree:::as.tree(tree, rf_output_total)

#####################
##### 2. Apply random forest to find the variants selected by VSURF using Rej vs. NoRej
#####################
load("/Users/Pinedasans/Catalyst/Results/ResultsRF_RejvsNoRej.Rdata")

dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])
phenotype<-as.factor(ifelse(demographics$phenotype[non.list]=="AMR","Rej",
                            ifelse(demographics$phenotype[non.list]=="CMR","Rej","NoRej")))

set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-phenotype[-i]
  testDataClass<-phenotype[i]

  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData,type="prob")[1,2]) # Prediction
}
library(ROCR)
predictions=as.vector(result_predicted)
pred=prediction(predictions,phenotype)

perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))


rf_output_total<- randomForest(phenotype~.,data=dataset,proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, phenotype)
legend("topright", legend=c("AMR+CMR","NoRej"),fill=brewer.pal(2,"Set1")) 


exome_variants_rf_selected_rej<-exome_variants_diff_complete[,fit$varselect.interp]
Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_rf_selected_rej)){
  Diff.AMR[i]<-table(exome_variants_rf_selected_rej[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_rf_selected_rej[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_rf_selected_rej[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_rf_selected_rej),df_joint_qc$snp_id)
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,rf_output_total$importance)
write.table(df_joint_sign,file="/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF_RejvsNoRej.txt",sep="\t",row.names = F)

library("pheatmap")
race <- ifelse(demographics$mismatch[non.list] == "0","Non race-mismatch", "race-mismatch")
related <- ifelse(demographics$REL2[non.list] == "0","Non related", "related")

annotation_row = data.frame(
  phenotype = phenotype,
  race = race,
  related = related)
rownames(annotation_row)<-variant_list$Id2
rownames(exome_variants_rf_selected_rej)<-variant_list$Id2
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

ann_colors = list (phenotype = c(Rej = "goldenrod","NoRej" = "darkolivegreen4"),
                   related = c("Non related" = "yellow", related = "lightblue"))
pheatmap(exome_variants_rf_selected_rej,cluster_rows = T,color = colorRampPalette(c("white", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors)


############################################
#### Apply VSURF in a two step process #####
###########################################

p.value<-read.table(file="/Users/Pinedasans/Catalyst/Results/p.value.endpoints.txt")
p.value<-p.value[,1]
exome_variants_sign<-exome_variants_diff2[,which(p.value<0.05)] #8,182

##Prepare data to run VSURF in the server
exome_variants_diff_complete<-exome_variants_sign[ , ! apply( exome_variants_sign, 2 , function(x) any(is.na(x)) ) ] #7,682
dataset<-data.frame(exome_variants_diff_complete)
save(dataset,demographics,file="/Users/Pinedasans/Data/Catalyst/DataRF_preselected.Rdata")

####Read the results after running in the server
load("/Users/Pinedasans/Catalyst/Results/ResultsRF_preselected.Rdata")
exome_variants_rf_selected_pre<-exome_variants_diff_complete[,fit$varselect.interp]
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected_pre),proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, demographics$phenotype[non.list])
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),
       fill=brewer.pal(length(levels(demographics$phenotype[non.list])),"Set1"))

# Create model with the variants found with the VSURF
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",3,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,1))
dataset<-data.frame(exome_variants_rf_selected_pre)

set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-demographics$phenotype[non.list][-i]
  testDataClass<-demographics$phenotype[non.list][i]
  
  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData)) # Prediction
}
result_roc <- multiclass.roc(endpoint, result_predicted)$auc


###Apply logistic regression to find only those that are risk 
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",1,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,3))
OR<-NULL
p.value.OR<-NULL
for (i in 1:ncol(exome_variants_rf_selected_pre)){
  print(i)
  model<-glm(exome_variants_rf_selected_pre[,i]~endpoint,family = "binomial")
  #m<-polr(factor(endpoint) ~ exome.variants.sign[,i],Hess=TRUE)
  OR[i]<-exp(coef(model))[2]
  p.value.OR[i]<-coefficients(summary(model))[2,4]
}

Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_rf_selected_pre)){
  Diff.AMR[i]<-table(exome_variants_rf_selected_pre[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_rf_selected_pre[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_rf_selected_pre[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_rf_selected_pre),df_joint_qc$snp_id)
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,OR,p.value.OR,rf_output$importance)
write.table(df_joint_sign,file="/Users/Pinedasans/Catalyst/Results/ResultsRF_2steps.txt",sep="\t",row.names = F)

library("pheatmap")
race <- ifelse(demographics$mismatch[non.list] == "0","Non race-mismatch", "race-mismatch")
related <- ifelse(demographics$REL2[non.list] == "0","Non related", "related")

annotation_row = data.frame(
  phenotype = demographics$phenotype[non.list],
  race = race,
  related = related)
rownames(annotation_row)<-variant_list$Id2
rownames(exome_variants_rf_selected_pre)<-variant_list$Id2
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

ann_colors = list (phenotype = c(AMR = "goldenrod", CMR = "darkorange", "No-REJ" = "darkolivegreen4"),
                   related = c("Non related" = "yellow", related = "lightblue"))
pheatmap(exome_variants_rf_selected_pre,cluster_rows = T,color = colorRampPalette(c("white", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors)


############################
#### Permutation Test RF ##
##########################
load("/Users/Pinedasans/Catalyst/Results/ResultsRF_permutation.Rdata")
p.value<-read.table(file="/Users/Pinedasans/Catalyst/Results/p.value.endpoints.txt")
p.value<-p.value[,1]
exome_variants_sign<-exome_variants_diff2[,which(p.value<0.05)] #8,182
exome_variants_diff_complete<-exome_variants_sign[ , ! apply( exome_variants_sign, 2 , function(x) any(is.na(x)) ) ] #7,682

var_select<-NULL
OOB<-NULL
for (i in 1:10){
  fit<-get(paste("fit",i,sep=""))
  exome_variants_rf_selected_perm<-exome_variants_diff_complete[,fit$varselect.inter]
  var_select[i]<-dim(exome_variants_rf_selected_perm)[2]
  rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected_perm),proximity=TRUE, keep.forest=T,ntree=100)
  OOB[i]<-rf_output_total$err.rate[100,1]
  tiff(paste("/Users/Pinedasans/Catalyst/Article/MDSplot_perm.tiff",i,".tiff",sep=""), width = 6, height =6, units = 'in', res = 300, compression = 'lzw')
    COLOR=brewer.pal(3,"Set1")
    MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = paste("Permutation",i,sep=""),cex=1.6)
    legend("bottomleft", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20,pt.cex=1.6 ,cex=1.2)
  dev.off()
}
    


########################
##### HLA analysis ####
#######################
id.hla.A<-grep("^HLA-A$",df_joint_qc$Gene.refGene)
id.hla.B<-grep("^HLA-B$",df_joint_qc$Gene.refGene)
id.hla.C<-grep("^HLA-C$",df_joint_qc$Gene.refGene)
id.hla.DPB1<-grep("^HLA-DPB1$",df_joint_qc$Gene.refGene)
id.hla.DQB1<-grep("^HLA-DQB1$",df_joint_qc$Gene.refGene)
id.hla.DRB1<-grep("^HLA-DRB1$",df_joint_qc$Gene.refGene)
id.hla.DRB5<-grep("^HLA-DRB5$",df_joint_qc$Gene.refGene)

id.HLA<-grep("HLA",df_joint_qc$Gene.refGene)
df_joint_hla<-df_joint_qc[id.HLA,]
id.hla.exome<-match(df_joint_hla$snp_id,colnames(exome_variants_diff_complete))
exome_variants_diff_hla<-exome_variants_diff_complete[,na.omit(id.hla.exome)]
save(exome_variants_diff_hla,demographics,file="/Users/Pinedasans/Catalyst/Data/DataRF_HLA.Rdata")

load("/Users/Pinedasans/Catalyst/Results/ResultsRF_HLA.Rdata")

dataset<-data.frame(exome_variants_diff_hla[,fit$varselect.interp])
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",3,
                 ifelse(demographics$phenotype[non.list]=="CMR",2,1))
set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-demographics$phenotype[non.list][-i]
  testDataClass<-demographics$phenotype[non.list][i]
  
  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData)) # Prediction
}
result_roc <- multiclass.roc(endpoint, result_predicted)$auc

exome_variants_rf_hla<-exome_variants_diff_hla[,fit$varselect.interp]
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_hla),proximity=TRUE, keep.forest=T,ntree=100)
png("/Users/Pinedasans/Catalyst/Results/MDSplot_hla.png")
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "HLA-region")
legend("topleft", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20)
dev.off()

########################################
### Association with Race mismatch  ####
########################################
load("/Users/Pinedasans/Catalyst/Results/ResultsRF_race.Rdata")

# Create model with the variants found with the VSURF
dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])
phenotype<-as.factor(ifelse(demographics$phenotype[non.list]=="AMR","Rej",
                            ifelse(demographics$phenotype[non.list]=="CMR","Rej","NoRej")))

set.seed(2)
result_roc<-NULL
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-phenotype[-i]
  testDataClass<-phenotype[i]
  
  rf_output <- randomForest(trainDataClass~.,data=trainData,proximity=TRUE, keep.forest=T,ntree=100)
  result_predicted <- c(result_predicted,predict(rf_output, testData,type="prob")[1,2]) # Prediction
}
library(ROCR)
predictions=as.vector(result_predicted)
pred=prediction(predictions,phenotype)

perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]

perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))


rf_output_total<- randomForest(phenotype~.,data=dataset,proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, phenotype)
legend("topright", legend=c("AMR+CMR","NoRej"),fill=brewer.pal(2,"Set1")) 




####Apply RF with the 19 and 13 variants overlapping
resultFisher<-read.table("/Users/Pinedasans/Catalyst/Data/Genotyping/13variantsOverlapping.txt",header=T,sep="\t")
id.snp<-match(resultFisher$SNP_rs,df_joint_qc$snp138)
id.fit<-match(df_joint_qc[id.snp,14],colnames(exome_variants_diff_complete))
exome_variants_rf_selected<-exome_variants_diff_complete[,na.omit(id.fit)] 
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected),proximity=TRUE, keep.forest=T,ntree=100)

##plot MDS
#png("/Users/Pinedasans/Catalyst/Results/MDSplot_orig.png")
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "13 overlap variants")
legend("bottomright", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20)
#dev.off()

####Apply RF with the variants selected by FisherExactTest
resultFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTest.txt",header=T,sep="\t")
id.fit<-match(resultFisher$snp_id,colnames(exome_variants_diff_complete))
exome_variants_rf_selected<-exome_variants_diff_complete[,na.omit(id.fit)] 
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=data.frame(exome_variants_rf_selected),proximity=TRUE, keep.forest=T,ntree=100)

##plot MDS
#png("/Users/Pinedasans/Catalyst/Results/MDSplot_orig.png")
COLOR=brewer.pal(3,"Set1")
MDSplot(rf_output_total, demographics$phenotype[non.list],palette = COLOR,main = "19 overlap variants")
legend("bottomright", legend=levels(demographics$phenotype[non.list]),col=COLOR, pch = 20)
#dev.off()

####################################
####Heatmap with the results#######
##################################
library("pheatmap")
race <- ifelse(demographics$mismatch[non.list] == "0","Non race-mismatch", "race-mismatch")
related <- ifelse(demographics$REL2[non.list] == "0","Non related", "related")

annotation_row = data.frame(
  phenotype = demographics$phenotype[non.list],
  race = race,
  related = related)
rownames(annotation_row)<-variant_list$Id2
rownames(exome_variants_rf_selected)<-variant_list$Id2
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
fill=brewer.pal(3,"Set1")
ann_colors = list (phenotype = c(AMR = fill[1], CMR = fill[2], "No-REJ" = fill[3]),
                   related = c("Non related" = "yellow", related = "lightblue"))
pheatmap(exome_variants_rf_selected ,cluster_rows = T,color = colorRampPalette(c("white", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors,border_color=F)
