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

#Work Directory
setwd("/Users/Pinedasans/Documents/Catalyst/Results/")
load("/Users/Pinedasans/Data/Catalyst/ExomeSeq_Diff_demo.Rdata")

#########################
##### Random Forest #### 
########################

##Prepare data to run VSURF in the server
exome_variants_diff_complete<-exome_variants_diff2[ , ! apply( exome_variants_diff2, 2 , function(x) any(is.na(x)) ) ] #450,981
save(exome_variants_diff_complete,demographics,file="/Users/Pinedasans/Data/Catalyst/DataRF.Rdata")


### 1. Apply random forest to find the variants selected by VSURF using three categories (AMR, CMR, NoRej)
load("/Users/Pinedasans/Documents/Catalyst/Results/ResultsRF.Rdata")

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

rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=dataset,proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, demographics$phenotype[non.list])
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),
       fill=brewer.pal(length(levels(demographics$phenotype[non.list])),"Set1")) 


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
write.table(df_joint_sign,file="/Users/Pinedasans/Documents/Catalyst/Results/ResultsEndpointRF.txt",sep="\t",row.names = F)

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

ann_colors = list (phenotype = c(AMR = "goldenrod", CMR = "darkorange", "No-REJ" = "darkolivegreen4"),
                   related = c("Non related" = "yellow", related = "lightblue"))
pheatmap(exome_variants_rf_selected,cluster_rows = T,color = colorRampPalette(c("green", "red"))(50),
         annotation_row = annotation_row, clustering_callback = callback,annotation_colors = ann_colors)

##### 2. Apply random forest to find the variants selected by VSURF using Rej vs. NoRej
load("/Users/Pinedasans/Documents/Catalyst/Results/ResultsRF_RejvsNoRej.Rdata")

dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])
phenotype<-as.factor(ifelse(demographics$phenotype[non.list]=="AMR","Rej",
                            ifelse(demographics$phenotype[non.list]=="CMR","Rej","NoRej")))

set.seed(2)
result_roc<-NULL
dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])
phenotype<-as.factor(ifelse(demographics$phenotype[non.list]=="AMR","Rej",
                            ifelse(demographics$phenotype[non.list]=="CMR","Rej","NoRej")))
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-phenotype[-i]
  testDataClass<-phenotype[i]
  
  #fit<-VSURF(x = trainData, trainDataClass,parallel = TRUE,ncores=2)
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


exome_variants_rf_selected<-exome_variants_diff_complete[,fit$varselect.interp]
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
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,rf_output$importance)
write.table(df_joint_sign,file="/Users/Pinedasans/Documents/Catalyst/Results/ResultsEndpointRF_RejvsNoRej.txt",sep="\t",row.names = F)


### 3. Apply random forest to find the variants selected by VSURF using AMR vs. CMRandNoRej
load("/Users/Pinedasans/Documents/Catalyst/Results/ResultsRF_AMRvsCMRNoRej.Rdata")

dataset<-data.frame(exome_variants_diff_complete[,fit$varselect.interp])
phenotype<-as.factor(ifelse(demographics$phenotype[non.list]=="CMR","NoRej","AMR"))

set.seed(2)
result_predicted<-NULL
for (i in 1:28){
  print(i)
  trainData <- dataset[-i,]
  testData <- dataset[i,]
  trainDataClass<-phenotype[-i]
  testDataClass<-phenotype[i]
  
  #fit<-VSURF(x = trainData, trainDataClass,parallel = TRUE,ncores=2)
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

rf_output_total <- randomForest(phenotype~.,data=dataset,proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, phenotype)
legend("topright", legend=c("AMR","CMR+NoRej"),fill=brewer.pal(2,"Set1")) 


############################################
#### Apply VSURF in a two step process #####
###########################################

p.value<-read.table(file="/Users/Pinedasans/Documents/Catalyst/Results/p.value.endpoints.txt")
exome_variants_sign<-exome_variants_diff2[,which(p.value<0.05)] #8,182

##Prepare data to run VSURF in the server
exome_variants_diff_complete<-exome_variants_sign[ , ! apply( exome_variants_sign, 2 , function(x) any(is.na(x)) ) ] #7,682
dataset<-data.frame(exome_variants_diff_complete)
save(dataset,demographics,file="/Users/Pinedasans/Data/Catalyst/DataRF_preselected.Rdata")

####Read the results after running in the server
load("/Users/Pinedasans/Documents/Catalyst/Results/ResultsRF_preselected.Rdata")
rf_output_total <- randomForest(demographics$phenotype[non.list]~.,data=dataset[,fit$varselect.interp],proximity=TRUE, keep.forest=T,ntree=100)
MDSplot(rf_output_total, demographics$phenotype[non.list])
legend("bottomleft", legend=levels(demographics$phenotype[non.list]),
       fill=brewer.pal(length(levels(demographics$phenotype[non.list])),"Set1"))
