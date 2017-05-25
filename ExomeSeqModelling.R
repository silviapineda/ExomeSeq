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
library("RColorBrewer")

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

####Functioanl variants
id<-match(rownames(exome_variants_present),df_joint_qc$snp_id)
df_joint_qc_diff<-df_joint_qc[id,]
table(df_joint_qc_diff$ExonicFunc.refGene)
counts = (matrix(data = c(19, 123, 59718, 474235), nrow = 2))
chisq.test(counts) #p-value = 0.5

df_joint_qc_diff_exonic<-df_joint_qc_diff[which(df_joint_qc_diff$ExonicFunc.refGene=="nonsynonymous SNV"),]


###########################################################################################
#### Obtain the main plot considering all the differences by pair ordered by endpoint  ####
###########################################################################################
y <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)

non.list<- seq(1,56,2) ##Donors

m <- rbind(c(2,3,4),c(1,1,1))
layout(m)
layout.show(4)

fill=brewer.pal(3,"Set1")

tiff("/Users/Pinedasans/Catalyst/Article/Boxplot.tiff", width = 6, height = 4, units = 'in', res = 300, compression = 'lzw')
boxplot(variant_mismatch~variant_list$phenotype,frame.plot = FALSE,col=fill,ylab="Variants Mismatched",
        ylim=c(40000,140000),cex=1.6)
dev.off()

num.AMR<-variant_mismatch[order(variant_list$phenotype)][1:14]
tiff("/Users/Pinedasans/Catalyst/Article/plotAMR.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'lzw')
plot(num.AMR[order(num.AMR,decreasing = T)],type="h", col=fill[1],ylim=c(40000,140000),ylab="Variants Mismatched",lty=1,lwd=10,xaxt="n",xlab="Nº pairs",xlim=c(0,15),
     main="AMR (386,958 mismatched variants)",frame.plot = FALSE)
axis(1,at= 1:14,label= paste("pair",rownames(variant_list),sep="")[order(variant_list$phenotype)][1:14][order(num.AMR,decreasing = T)],las=2,cex.axis=0.9)
text(num.AMR[order(num.AMR,decreasing = T)]+4000,labels=demographics$RACE[y==1][order(variant_list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.8)
text(num.AMR[order(num.AMR,decreasing = T)]+8000,labels=demographics$RACE[y==0][order(variant_list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.8)
text(num.AMR[order(num.AMR,decreasing = T)]+12000,labels=demographics$RELplot[y==0][order(variant_list$phenotype)][1:14][order(num.AMR,decreasing = T)],cex=0.8)
dev.off()


num.CMR<-variant_mismatch[order(variant_list$phenotype)][15:21]
tiff("/Users/Pinedasans/Catalyst/Article/plotCMR.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'lzw')
plot(num.CMR[order(num.CMR,decreasing = T)],type="h", col=fill[2],ylim=c(40000,140000),ylab="Variants Mismatched",lty=1,lwd=10,xaxt="n",xlab="Nº pairs",xlim=c(0,8),
     main="CMR (268,722 mismatched variants)",frame.plot = FALSE)
axis(1,at= 1:7,label= paste("pair",rownames(variant_list),sep="")[order(variant_list$phenotype)][15:21][order(num.CMR,decreasing = T)],las=2,cex.axis=0.9)
text(num.CMR[order(num.CMR,decreasing = T)]+4000,labels=demographics$RACE[y==1][order(variant_list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.8)
text(num.CMR[order(num.CMR,decreasing = T)]+8000,labels=demographics$RACE[y==0][order(variant_list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.8)
text(num.CMR[order(num.CMR,decreasing = T)]+12000,labels=demographics$RELplot[y==0][order(variant_list$phenotype)][15:21][order(num.CMR,decreasing = T)],cex=0.8)
dev.off()

num.NoRej<-variant_mismatch[order(variant_list$phenotype)][22:28]
tiff("/Users/Pinedasans/Catalyst/Article/plotNoRej.tiff", width = 5, height = 5, units = 'in', res = 300, compression = 'lzw')
plot(num.NoRej[order(num.NoRej,decreasing = T)],type="h", col=fill[3],ylim=c(40000,140000),ylab="Variants Mismatched",lwd=10,lty=1,xaxt="n",xlab="Nº pairs",xlim=c(0,8),
     main="NoRej (248,531 mismatched variants)",frame.plot = FALSE)
axis(1,at= 1:7,label= paste("pair",rownames(variant_list),sep="")[order(variant_list$phenotype)][22:28][order(num.NoRej,decreasing = T)],las=2,cex.axis=0.9)
text(num.NoRej[order(num.NoRej,decreasing = T)]+4000,labels=demographics$RACE[y==1][order(variant_list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.8)
text(num.NoRej[order(num.NoRej,decreasing = T)]+8000,labels=demographics$RACE[y==0][order(variant_list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.8)
text(num.NoRej[order(num.NoRej,decreasing = T)]+12000,labels=demographics$RELplot[y==0][order(variant_list$phenotype)][22:28][order(num.NoRej,decreasing = T)],cex=0.8)
dev.off()

summary(lm(variant_mismatch~demographics$phenotype[non.list]+distance[,1]))

distance<-read.table("/Users/Pinedasans/Catalyst/Results/Distance.PCA.txt")
summary(lm(distance[,1]~demographics$phenotype[non.list]))
boxplot(distance[,1]~variant_list$phenotype,frame.plot = FALSE,col=fill,ylab="Variants Mismatched",
      cex=1.2)
plot(distance[,1],variant_mismatch)

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

##45% of the times we have significance


####################################################################
########## MODELLING THE DIFFERENCES ###############################
####################################################################

##Looking for significance between the variant differed and endpoint
endpoint<-ifelse(demographics$phenotype[non.list]=="No-REJ",0,
                 ifelse(demographics$phenotype[non.list]=="CMR",1,2))

##Preparing the data with the diff
# exome_variants_present_t<-t(exome_variants_present)
# exome_variants_diff<-apply(exome_variants_present_t,2,function(x) abs(diff(x))[non.list])
# exome_variants_diff2<-apply(exome_variants_diff,2,function(x) replace(x,x==2,1))
# save(exome_variants_diff2,demographics,variant_list,file="/Users/Pinedasans/Catalyst/Data/ExomeSeq_Diff_demo.Rdata")

load("/Users/Pinedasans/Catalyst/Data/ExomeSeq_Diff_demo.Rdata")

###Considering a fisher.test
p.value<-rep(NA,ncol(exome_variants_diff2))
for (i in 1:ncol(exome_variants_diff2)){
  print(i)
  tab<-table(exome_variants_diff2[,i],demographics$phenotype[non.list])
  if(dim(tab)[1]>1){
    p.value[i]<-fisher.test(tab)$p.value
  }
}

write.table(p.value,file="/Users/Pinedasans/Catalyst/Results/p.value.endpoints.txt")
p.value<-read.table(file="/Users/Pinedasans/Catalyst/Results/p.value.endpoints.txt")
p.value<-p.value[,1]
p.value.adj<-p.adjust(p.value,method = "BH") #There are no significant results after MT correction
exome_variants_sign<-exome_variants_diff2[,which(p.value<0.001)] #123
p.value.sign<-p.value[which(p.value<0.001)]

##Plotting manhattan plot
library("qqman")
id.man<-match(colnames(exome_variants_diff2),df_joint_qc$snp_id)
data.manhattan<-cbind(df_joint_qc[na.omit(id.man),c(14,1,2)],p.value)
colnames(data.manhattan)=c("SNP","CHR","BP","P")
data.manhattan$CHR<-as.numeric(as.character(data.manhattan$CHR))
data.manhattan.NoNan<-data.manhattan[which(is.na(data.manhattan$P)==F),]
data.manhattan.order<-data.manhattan.NoNan[order(data.manhattan.NoNan$CHR,data.manhattan.NoNan$BP),]
highlightSNPs <- data.manhattan.order[which(data.manhattan.order$P<=0.001),1]
grep("1:26671084_G/T",data.manhattan.order$SNP)
snpsOfInterest<-data.manhattan.order$SNP[11100:11800]
tiff("/Users/Pinedasans/Catalyst/Article/Manhattan.tiff", width = 14, height = 8, units = 'in', res = 300, compression = 'lzw')
manhattan(data.manhattan.NoNan, ylim = c(0, 6),highlight=highlightSNPs,suggestiveline = FALSE,genomewideline = FALSE,
          cex=1.2)
dev.off()
#abline(h = -log10(0.001), col = "blue")
#data locus
dataLocus<-cbind(df_joint_qc[na.omit(id.man),c(1,2,2,14)],p.value)
colnames(dataLocus)=c("CHROM",	"BEGIN",	"END",	"MARKER_ID","PVALUE")
write.table(dataLocus,"/Users/Pinedasans/Catalyst/Results/dataLocus.txt",row.names = F,sep="\t")

####
#make annotation factor
ann<-rep(1, length(data.manhattan.NoNan$P))
ann[with(data.manhattan.NoNan, CHR==1 & BP>= 227057885 & BP<227083806)]<-2
ann[with(data.manhattan.NoNan, CHR==8 & BP>=113235157 & BP<114449328)]<-3
ann[with(data.manhattan.NoNan, CHR==19 & BP>=2100988 & BP<2164464)]<-4
ann[with(data.manhattan.NoNan, CHR==10 & BP>=12237964 & BP<12292588)]<-5
ann[with(data.manhattan.NoNan, CHR==1 & BP>=26648350 & BP<26680621)]<-6
ann[with(data.manhattan.NoNan, CHR==1 & BP>=26648350 & BP<26680621)]<-7
ann[with(data.manhattan.NoNan, CHR==20 & BP>=36838890 & BP<36889174)]<-8
ann[with(data.manhattan.NoNan, CHR==9 & BP>=130267618 & BP<130341268)]<-9
ann[with(data.manhattan.NoNan, CHR==7 & BP>=100547257 & BP<100550424)]<-10
ann[with(data.manhattan.NoNan, CHR==11 & BP>=4790209 & BP<4791168)]<-11
ann[with(data.manhattan.NoNan, CHR==11 & BP>=124120423 & BP<124135763)]<-12
ann[with(data.manhattan.NoNan, CHR==11 & BP>=124134723 & BP<124135763)]<-13

ann<-factor(ann, levels=1:13, labels=c("","PSEN2","CSMD3","AP3D1","CDC123","AIM1L","CHRNA10","KIAA1755",
                                      "FAM129B","MUC3A","OR51F1","OR8G1","OR8G5"))
#draw plot with annotation
manhattan.plot(data.manhattan.NoNan$CHR,data.manhattan.NoNan$BP, data.manhattan.NoNan$P,
               annotate=list( ann, "PSEN2"=list(col=fill[1]),"CSMD3"=list(col=fill[1]),
                              "AP3D1"=list(col=fill[1]),"CDC123"=list(col=fill[1]),
                              "AIM1L"=list(col=fill[2]),"CHRNA10"=list(col=fill[2]),
                              "KIAA1755"=list(col=fill[2]),"FAM129B"=list(col=fill[1]),
                              "MUC3A"=list(col=fill[1]),"OR51F1"=list(col=fill[1]),
                              "OR8G1"=list(col=fill[1]),"OR8G5"=list(col=fill[1])))

png("mymanhattan.png", width=950, height=500)
print(manhattan.plot())
dev.off()


#############################
### Permutation Analysis ###
############################
###Considering a fisher.test
p.value<-matrix(NA,100,ncol(exome_variants_diff2))
for (b in 1:100){
  print(b)
  perm<-sample(demographics$phenotype[non.list])
  for (i in 1:ncol(exome_variants_diff2)){
    tab<-table(exome_variants_diff2[,i],perm)
    if(dim(tab)[1]>1){
      p.value[b,i]<-fisher.test(tab)$p.value
    }
  }
}
write.csv(p.value,"/Users/Pinedasans/Catalyst/Results/p_value_perm_fisher.csv")
colnames(p.value)<-colnames(exome_variants_diff2)
numVar<-NULL
for (i in 1:100){
  numVar[i]<-table(p.value[i,]<0.001)[2]
}
p.value.sign.perm<-p.value[,match(colnames(exome_variants_sign),colnames(p.value))]
numSign<-NULL
for (i in 1:ncol(p.value.sign.perm)){
  numSign[i]<-table(p.value.sign.perm[,i]<0.05)[2]
}

library(MASS)
#endpoint <- relevel(factor(endpoint), ref = "0")
OR<-NULL
p.value.OR<-NULL
for (i in 1:ncol(exome_variants_sign)){
  print(i)
  #model<-glm(exome_variants_sign[,i]~endpoint,family = "binomial")
  m<-polr(factor(endpoint) ~ exome_variants_sign[,i],Hess=TRUE)
  OR[i]<-exp(coef(summary(m)))[1,1]
  p.value.OR[i]<- pnorm(abs(coef(summary(m))[1, "t value"]), lower.tail = FALSE) * 2
 #p.value.OR[i]<-coefficients(summary(model))[2,4]
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

###All differences to plot in the circos plot
Diff.AMR<-NULL
Diff.CMR<-NULL
Diff.NoRej<-NULL
for(i in 1:ncol(exome_variants_diff2)){
  Diff.AMR[i]<-table(exome_variants_diff2[which(demographics$phenotype[non.list]=="AMR"),i]!=0)['TRUE']
  Diff.CMR[i]<-table(exome_variants_diff2[which(demographics$phenotype[non.list]=="CMR"),i]!=0)['TRUE']
  Diff.NoRej[i]<-table(exome_variants_diff2[which(demographics$phenotype[non.list]=="No-REJ"),i]!=0)['TRUE']
}
###Annotate the variants
id.joint<-match(colnames(exome_variants_diff2),df_joint_qc$snp_id)
df_joint_all<-cbind(df_joint_qc[id.joint,],Diff.AMR,Diff.CMR,Diff.NoRej)
df_joint_all<-df_joint_all[,c(1,2,3,72,73,74)]

####To plot in the circos plot
for (i in 1:22){
  df_joint_all$Chr<-replace(as.character(df_joint_all$Chr),as.character(df_joint_all$Chr)==i,paste("hs",i,sep=""))
}
df_joint_all$Diff.AMR<-replace(df_joint_all$Diff.AMR,is.na(df_joint_all$Diff.AMR)==T,0)
df_joint_all$Diff.CMR<-replace(df_joint_all$Diff.CMR,is.na(df_joint_all$Diff.CMR)==T,0)
df_joint_all$Diff.NoRej<-replace(df_joint_all$Diff.NoRej,is.na(df_joint_all$Diff.NoRej)==T,0)
write.table(df_joint_all,"/Users/Pinedasans/programs/circos-0.69-5/catalyst/variants.hist.all.txt")


###Count the variants in the donor and the variants in the recipient
df_joint_sign[,16:71]

donor_res<-rep(NA,ncol(exome_variants_sign))
recipient_res<-rep(NA,ncol(exome_variants_sign))
match_res<-rep(NA,ncol(exome_variants_sign))
for (i in 1:ncol(exome_variants_sign)){
  j=16
  donor<-0
  recipient<-0
  match<-0
  while(j<72){
    print(j)
    if(df_joint_sign[i,j]!="NA" & df_joint_sign[i,(j+1)]!="NA"){
      if(as.numeric(df_joint_sign[i,j])-as.numeric(df_joint_sign[i,(j+1)])>0){
        donor=donor+1
        j=j+2
      } else if(as.numeric(df_joint_sign[i,j])-as.numeric(df_joint_sign[i,(j+1)])<0){
          recipient=recipient+1
          j=j+2
      } else{
        match=match+1
        j=j+2
      }
    } else {
      j=j+2
    }
  }
  donor_res[i]<-donor
  recipient_res[i]<-recipient
  match_res[i]<-match
}

write.table(cbind(df_joint_sign[,-c(16:71)],donor_res,recipient_res,match_res),file="/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTestSign.txt",sep="\t",row.names = F)




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
id.hla.DQA1<-grep("^HLA-DQA1$",df_joint_qc$Gene.refGene)

id.HLA<-grep("HLA",df_joint_qc$Gene.refGene)

####Analysis 1: Association analysis for each SNP in HLA region and clinical endpoint

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
p.value.HLA.DQA1<-FisherTestSNV(id.hla.DQA1,df_joint_qc)

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

####Applying multinomial adjusted by race
library(MASS)
OR<-NULL
OR_race<-NULL
p.value.OR<-NULL
p.value.OR.race<-NULL
for (i in 1:ncol(exome_variants_diff2)){
  print(i)
  m<-polr(factor(endpoint) ~ exome_variants_diff2[,i] + demographics$mismatch[non.list],Hess=TRUE)
  OR[i]<-exp(coef(summary(m)))[1,1]
  OR_race[i]<-exp(coef(summary(m)))[1,2]
  p.value.OR[i]<- pnorm(abs(coef(summary(m))[1, "t value"]), lower.tail = FALSE) * 2
  p.value.OR.race[i]<- pnorm(abs(coef(summary(m))[2, "t value"]), lower.tail = FALSE) * 2
}
write.table(cbind(OR,p.value.OR,OR_race,p.value.OR.race),file="/Users/Pinedasans/Catalyst/Results/ResultsMultinomialAdjByRace.txt")

p.value.OR.2<-p.value.OR[which(is.na(p.value.OR)==FALSE)]
p.value.OR.3<-p.value.OR.2[which(p.value.OR.2!=0)]


###########################################################################################################################################################
