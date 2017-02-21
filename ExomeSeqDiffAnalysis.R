
rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Esome Sequencing Donor/Recipient pairs
###
### CITATION: 
###
### PROCESS: Analyzing ExomeSeq Data
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
setwd("/Users/Pinedasans/Data/Catalyst/ExomeSeq/VCF/")


###Reading the Total genotypes data frame
load("ExomeSeqVCF.Rdata")

###Read the non-annotated variants 
# df.non.annotated<-read.table("/Users/Pinedasans/Data/Catalyst/ExomeSeq/annovar/non.annotated.variants.txt")
# id.non.annotated<-match(df.non.annotated$snp_id,df.joint$snp_id)
# df.joint<-df.joint[-id.non.annotated,] ### Total variants to analyze 513,341 

##To trasnform into numeric matrix
exome.variants <- apply(df.joint[,6:60],2,as.numeric)
rownames(exome.variants) <- df.joint$snp_id

##QC of missingness
exome.variants.qc<-exome.variants[rowSums(is.na(exome.variants) == FALSE) >= 52,] #486,089 that has at least 95% of non-missingness
id.qc<-match(rownames(exome.variants.qc),df.joint$snp_id)
df.joint.qc<-df.joint[na.omit(id.qc),]


save(df.joint.qc,exome.variants.qc,file="ExomeSeqVCF_QC.Rdata")

####Use only those that are annotated
my.annovar <- read.csv("/Users/Pinedasans/Data/Catalyst/ExomeSeq/annovar/myanno.joint.output.snp.indel.recal.QC.avinput.hg19_multianno_columnfilter.csv")
df.joint.annotated<-df.joint.qc
colnames(df.joint.annotated)[2:3]<-c("Start","Chr")
merge.data <- merge(my.annovar,df.joint.annotated, by=c("Chr","Start")) #486,203
merge.data <- merge.data[,c(1:68,69,71,70,71)]
colnames(merge.data)[c(70,72)]<-c("Discovery-270862-R1","Discovery-270862-R2")
id.variants.annotated<-match(df.joint.annotated$snp_id,merge.data$snp_id)
df.joint.annotated<-merge.data[na.omit(id.variants.annotated),] #485,984 annotated variants
exome.variants.annotated <- apply(df.joint.annotated[,17:72],2,as.numeric)
rownames(exome.variants.annotated) <- df.joint.annotated$snp_id

save(df.joint.annotated,exome.variants.annotated,file="ExomeSeqVCF_QC_annotated.Rdata")  


##Mirar a partir de aqui que pasa
##Looking for all the differences between pairs
j=1
for (i in 1:28){
  assign(paste("exome.pair",(i),sep=""), which(exome.variants.annotated[,j]!=exome.variants.annotated[,j+1]))
  assign(paste("df.joint.pair",(i),sep=""), df.joint.annotated[get(paste("exome.pair",(i),sep="")),c(1:16,j+16,j+17)])
  j=j+2
}

###Divide between the variants present only in the donor and the variants only present in  the recipient
for (i in 1:28){
  assign(paste("df.joint.pair",(i),".presentDonor",sep=""),
         get(paste("df.joint.pair",(i),sep=""))[which((get(paste("df.joint.pair",(i),sep=""))[,18]=="0" & get(paste("df.joint.pair",(i),sep=""))[,17]!="0") | 
                                                        (get(paste("df.joint.pair",(i),sep=""))[,18]=="1" & (get(paste("df.joint.pair",(i),sep=""))[,17]=="2"))),])
  assign(paste("df.joint.pair",(i),".presentRecipient",sep=""),
         get(paste("df.joint.pair",(i),sep=""))[which((get(paste("df.joint.pair",(i),sep=""))[,17]=="0" & get(paste("df.joint.pair",(i),sep=""))[,18]!="0") | 
                                                        (get(paste("df.joint.pair",(i),sep=""))[,17]=="1" & (get(paste("df.joint.pair",(i),sep=""))[,18]=="2"))),])
  
}


###Variants common in the donor and recipient
j=1
for (i in 1:28){
  assign(paste("exome.common.pair",(i),sep=""), which(exome.variants.annotated[,j]==exome.variants.annotated[,j+1]))
  assign(paste("df.joint.common.pair",(i),sep=""), df.joint.annotated[get(paste("exome.common.pair",(i),sep="")),c(1:16,j+16,j+17)])
  j=j+2
}

##Count number of variants per pair that are DIFFERENT
num.presentDonor <- NULL
num.presentRecipient <- NULL
for (i in 1:28) {
  num.presentDonor[i]<- dim(get(paste("df.joint.pair",(i),".presentDonor",sep="")))[1]
  num.presentRecipient[i]<-dim(get(paste("df.joint.pair",(i),".presentRecipient",sep="")))[1]
}


##Count number of variants per pair that are COMMON
num.common <- NULL
for (i in 1:28) {
  num.common[i]<- dim(get(paste("df.joint.common.pair",(i),sep="")))[1]
}

id<-seq(6,60,2)
xx<-cbind(num.presentDonor,num.presentRecipient,num.common)
rownames(xx)<-colnames(df.joint)[id]
write.table(xx,file="num.variants.annotated.txt",row.names = T,col.names = T,sep="\t")

##Save in an Rdata to work further
save(list = ls(all.names = TRUE), file = "DiffByPairsAnnotated.RData")

