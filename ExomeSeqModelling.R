rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Esome Sequencing Donor/Recipient pairs
###
### CITATION: 
###
### PROCESS: Analyzing the differences by Pair 
###          In this version we can consider the variants that have passed the QC threshold and are annotated
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

#Work Directory
setwd("/Users/Pinedasans/Catalyst/Results/")

###Reading the Total genotypes data frame
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq/ExomeSeqVCF_SNPs.Rdata")

##Load the data with the difference
load("/Users/Pinedasans/Catalyst/Data/ExomeSeq/DiffByPairs.RData")

###List the variants that are present in Donor and not in recipient by outcome
variant_list<-read.table("/Users/Pinedasans/Catalyst/Data/ExomeSeq/Total_variants.txt",header=T ) ##Reading the variant list with the phenotype

###List the Demographic info
demographics<-read.table("/Users/Pinedasans/Catalyst/Data/Demographics.txt",sep="\t",header=T)

#Demographics recipients
table(demographics$SEX[non.list],demographics$phenotype[non.list])
summary(demographics$Age[non.list][which(demographics$phenotype[non.list]=="AMR")]) 
summary(demographics$Age[non.list][which(demographics$phenotype[non.list]=="CMR")]) 
summary(demographics$Age[non.list][which(demographics$phenotype[non.list]=="No-REJ")]) 
sd(demographics$Age[non.list][which(demographics$phenotype[non.list]=="AMR")]) 
sd(demographics$Age[non.list][which(demographics$phenotype[non.list]=="CMR")]) 
sd(demographics$Age[non.list][which(demographics$phenotype[non.list]=="No-REJ")]) 
table(demographics$RACE[non.list],demographics$phenotype[non.list])
table(demographics$REL[non.list],demographics$phenotype[non.list])


##Count number of variants per pair that are DIFFERENT
variant_mismatch <- NULL
for (i in 1:28) {
  variant_mismatch[i]<- dim(get(paste("df_joint_pair",(i),sep="")))[1]
}

variantList_totalMismatch <- NULL
variantList_totalMismatch_AMR<-NULL
variantList_totalMismatch_CMR<-NULL
variantList_totalMismatch_NoRej<-NULL

for (i in 1:28){
  print(i)
  variantList_totalMismatch <- c(variantList_totalMismatch,as.character(get(paste("df_joint_pair",i,sep=""))[,14]))
}
for (i in order(variant_list$phenotype)[1:14]){
  print(i)
  variantList_totalMismatch_AMR <- c(variantList_totalMismatch_AMR,as.character(get(paste("df_joint_pair",i,sep=""))[,14]))
}
for (i in order(variant_list$phenotype)[15:21]){
  print(i)
  variantList_totalMismatch_CMR <- c(variantList_totalMismatch_CMR,as.character(get(paste("df_joint_pair",i,sep=""))[,14]))
}
for (i in order(variant_list$phenotype)[22:28]){
  print(i)
  variantList_totalMismatch_NoRej <- c(variantList_totalMismatch_NoRej,as.character(get(paste("df_joint_pair",i,sep=""))[,14]))
}

variantList_totalMismatch_unique <- unique(variantList_totalMismatch) #474,235
variantList_totalMismatch_AMR_unique <- unique(variantList_totalMismatch_AMR) #386,958
variantList_totalMismatch_CMR_unique <- unique(variantList_totalMismatch_CMR) #268,722
variantList_totalMismatch_NoRej_unique <- unique(variantList_totalMismatch_NoRej) #248,531


id.present <- match(variantList_totalMismatch_unique,rownames(exome_variants_qc))
id.present_AMR <- match(variantList_totalMismatch_AMR_unique,rownames(exome_variants_qc))
id.present_CMR <- match(variantList_totalMismatch_CMR_unique,rownames(exome_variants_qc))
id.present_NoRej<- match(variantList_totalMismatch_NoRej_unique,rownames(exome_variants_qc))


exome_variants_present <- exome_variants_qc[na.omit(id.present),]
exome_variants_present_AMR <- exome_variants_qc[na.omit(id.present_AMR),]
exome_variants_present_CMR <- exome_variants_qc[na.omit(id.present_CMR),]
exome_variants_present_NoRej <- exome_variants_qc[na.omit(id.present_NoRej),]


###########################################################################################
#### Obtain the main plot considering all the differences by pair ordered by endpoint  ####
###########################################################################################
y <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)

non.list<- seq(1,56,2) ##Donors

m <- rbind(c(2,3,4),c(1,1,1))
layout(m)
layout.show(4)

boxplot(variant_mismatch~variant_list$phenotype,frame.plot = FALSE,col=c("goldenrod","darkorange","darkolivegreen4"),ylab="Num Variants Differ",ylim=c(40000,140000))

num.AMR<-variant_mismatch[order(variant_list$phenotype)][1:14]
plot(num.AMR[order(num.AMR,decreasing = T)],type="h", col="goldenrod",ylim=c(40000,140000),ylab="Num Variants Differ",lty=1,lwd=5,xaxt="n",xlab="Nº pairs",
     main="AMR (386,958 variants)",frame.plot = FALSE)
axis(1,at= 1:14,label= paste("pair",rownames(variant_list),sep="")[order(variant.list$phenotype)][1:14][order(num.AMR,decreasing = T)],las=2,cex.axis=0.8)
text(num.AMR[order(num.AMR,decreasing = T)]+3000,labels=demographics$RACE[y==1][order(variant.list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.9)
text(num.AMR[order(num.AMR,decreasing = T)]+6000,labels=demographics$RACE[y==0][order(variant.list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.9)
text(num.AMR[order(num.AMR,decreasing = T)]+9000,labels=demographics$RELplot[y==0][order(variant.list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.9)

num.CMR<-variant_mismatch[order(variant_list$phenotype)][15:21]
plot(num.CMR[order(num.CMR,decreasing = T)],type="h", col="darkorange",ylim=c(40000,140000),ylab="Num Variants Differ",lty=1,lwd=5,xaxt="n",xlab="Nº pairs",
     main="CMR (268,722 variants)",frame.plot = FALSE)
axis(1,at= 1:7,label= paste("pair",rownames(variant_list),sep="")[order(variant.list$phenotype)][15:21][order(num.CMR,decreasing = T)],las=2,cex.axis=0.8)
text(num.CMR[order(num.CMR,decreasing = T)]+3000,labels=demographics$RACE[y==1][order(variant.list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.9)
text(num.CMR[order(num.CMR,decreasing = T)]+6000,labels=demographics$RACE[y==0][order(variant.list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.9)
text(num.CMR[order(num.CMR,decreasing = T)]+9000,labels=demographics$RELplot[y==0][order(variant.list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.9)

num.NoRej<-variant_mismatch[order(variant_list$phenotype)][22:28]
plot(num.NoRej[order(num.NoRej,decreasing = T)],type="h", col="darkolivegreen4",ylim=c(40000,140000),ylab="Num Variants Differ",lwd=5,lty=1,xaxt="n",xlab="Nº pairs",
     main="No-Rej (248,531 variants)",frame.plot = FALSE)
axis(1,at= 1:7,label= paste("pair",rownames(variant_list),sep="")[order(variant.list$phenotype)][22:28][order(num.NoRej,decreasing = T)],las=2,cex.axis=0.8)
text(num.NoRej[order(num.NoRej,decreasing = T)]+3000,labels=demographics$RACE[y==1][order(variant.list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.9)
text(num.NoRej[order(num.NoRej,decreasing = T)]+6000,labels=demographics$RACE[y==0][order(variant.list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.9)
text(num.NoRej[order(num.NoRej,decreasing = T)]+9000,labels=demographics$RELplot[y==0][order(variant.list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.9)


summary(lm(variant_mismatch~demographics$phenotype[non.list]))

##Selecting randomly 7 pairs from the AMR to find if there is significance with the same number of pairs in the AMR group
p.value<-NULL
for (b in 1:1000){
  xx<-sample(1:14,7,replace = T)
  yy<-sample(15:21,7,replace = T)
  zz<-sample(22:28,7,replace = T)
  model<-summary(lm(c(variant_mismatch[order(variant_list$phenotype)][xx],variant_mismatch[order(variant_list$phenotype)][yy],
                      variant_mismatch[order(variant_list$phenotype)][zz])~factor(c(demographics$phenotype[non.list][order(variant_list$phenotype)][xx],demographics$phenotype[non.list][order(variant_list$phenotype)][yy],
                                                                       demographics$phenotype[non.list][order(variant_list$phenotype)][zz]))))
  p.value[b]<-coef(model)[3,4]
}

##42% of the times we have significance


####################################################################
########## MODELLING THE DIFFERENCES ###############################
####################################################################

##Looking for significance between the variant differed and endpoint
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",0,
                 ifelse(demographics$phenotype[non.list]=="CMR",1,2))

exome_variants_present_t<-t(exome_variants_present)
exome_variants_diff<-apply(exome_variants_present_t,2,function(x) abs(diff(x))[non.list])
exome_variants_diff2<-apply(exome_variants_diff,2,function(x) replace(x,x==2,1))

save(exome_variants_diff2,demographics,variant_list,file="/Users/Pinedasans/Data/Catalyst/ExomeSeq_Diff_demo.Rdata")

###Considering a fisher.test
p.value<-rep(NA,ncol(exome_variants_diff2))
for (i in 1:ncol(exome_variants_diff2)){
  print(i)
  tab<-table(exome_variants_diff2[,i],demographics$phenotype[non.list])
  if(dim(tab)[1]>1){
    p.value[i]<-fisher.test(tab)$p.value
  }
}

write.table(p.value,file="/Users/Pinedasans/Documents/Catalyst/Results/p.value.endpoints.txt")
p.value<-read.table(file="/Users/Pinedasans/Documents/Catalyst/Results/p.value.endpoints.txt")
p.value<-p.value[,1]
p.value.adj<-p.adjust(p.value,method = "BH") #There are no significant results after MT correction
exome_variants_sign<-exome_variants_diff2[,which(p.value<0.001)] #123
p.value.sign<-p.value[which(p.value<0.001)]

OR<-NULL
p.value.OR<-NULL
for (i in 1:ncol(exome_variants_sign)){
  print(i)
  model<-glm(exome_variants_sign[,i]~endpoint,family = "binomial")
  #m<-polr(factor(endpoint) ~ exome.variants.sign[,i],Hess=TRUE)
  OR[i]<-exp(coef(model))[2]
  p.value.OR[i]<-coefficients(summary(model))[2,4]
}


#####
Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_sign)){
  Diff.AMR[i]<-table(exome_variants_sign[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_sign[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_sign[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_sign),df_joint_qc$snp_id)
df_joint_sign<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,p.value.sign,OR,p.value.OR)
write.table(df_joint_sign,file="/Users/Pinedasans/Documents/Catalyst/Results/ResultsEndpointFisherTestSign90%.txt",sep="\t",row.names = F)


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


FisherTestSNV<- function (gene.id,DATA){
###Considering a fisher.test
  DATA_HLA<-DATA[gene.id,]
  id.snv<-match(DATA_HLA$snp_id,colnames(exome_variants_diff2))
  exome_variants_diff_hla<-exome_variants_diff2[,na.omit(id.snv)]
  p.value.hla<-rep(NA,ncol(exome_variants_diff_hla))
  for (i in 1:ncol(exome_variants_diff_hla)){
    print(i)
    tab<-table(exome_variants_diff_hla[,i],demographics$phenotype[non.list])
    if(dim(tab)[1]>1){
      p.value.hla[i]<-fisher.test(tab)$p.value
    }
  }
  return(p.value.hla)
}

p.value.HLA<-FisherTestSNV(id.HLA,df_joint_qc)
p.value.HLA.A<-FisherTestSNV(id.hla.A,df_joint_qc)
p.value.HLA.B<-FisherTestSNV(id.hla.B,df_joint_qc)
p.value.HLA.C<-FisherTestSNV(id.hla.C,df_joint_qc)
p.value.HLA.DPB1<-FisherTestSNV(id.hla.DPB1,df_joint_qc)
p.value.HLA.DQB1<-FisherTestSNV(id.hla.DQB1,df_joint_qc)
p.value.HLA.DRB1<-FisherTestSNV(id.hla.DRB1,df_joint_qc)
p.value.HLA.DRB5<-FisherTestSNV(id.hla.DRB5,df_joint_qc)

exome_variants_HLA<-exome_variants_diff2[,id.HLA]

Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_HLA)){
  Diff.AMR[i]<-table(exome_variants_HLA[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_HLA[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_HLA[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}


###########
### Race mismatch
##########
non.list<-seq(1,56,2)
p.value.race<-rep(NA,ncol(exome_variants_diff2))
for (i in 1:ncol(exome_variants_diff2)){
  print(i)
  tab<-table(exome_variants_diff2[,i],demographics$mismatch[non.list])
  if(dim(tab)[1]>1){
    p.value.race[i]<-fisher.test(tab)$p.value
  }
}

write.table(p.value.race,file="/Users/Pinedasans/Documents/Catalyst/Results/p.value.race.txt")
p.value.race<-read.table(file="/Users/Pinedasans/Documents/Catalyst/Results/p.value.race.txt")
p.value.race<-p.value.race[,1]
p.value.race.adj<-p.adjust(p.value.race,method = "fdr")
exome_variants_race_sign<-exome_variants_diff2[,which(p.value.race.adj<0.1)] #1,206
p.value.race.sign<-p.value.race[which(p.value.race.adj<0.1)]

Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL

for(i in 1:ncol(exome_variants_race_sign)){
  Diff.AMR[i]<-table(exome_variants_race_sign[which(demographics$phenotype[non.list]=="AMR" & demographics$mismatch[non.list]==1),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_race_sign[which(demographics$phenotype[non.list]=="CMR" & demographics$mismatch[non.list]==1),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_race_sign[which(demographics$phenotype[non.list]=="No-REJ" & demographics$mismatch[non.list]==1),i]!=0)['TRUE']
}

###Annotate the variants
id.joint<-match(colnames(exome_variants_race_sign),df_joint_qc$snp_id)
df_joint_sign_race<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej,p.value.race.sign)
write.table(df_joint_sign_race,file="/Users/Pinedasans/Documents/Catalyst/Results/ResultsRaceFisherTest.txt")

###########################################################################################################################################################
