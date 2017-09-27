
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
setwd("/Users/Pinedasans/Catalyst/Data/ExomeSeq/")

###Reading the Total genotypes data frame
load("ExomeSeqVCF_SNPs.Rdata")

##Looking for all the differences between pairs
j=1
for (i in 1:28){
  assign(paste("exome_pair",(i),sep=""), which(exome_variants_qc[,j]!=exome_variants_qc[,j+1]))
  assign(paste("df_joint_pair",(i),sep=""), df_joint_qc[get(paste("exome_pair",(i),sep="")),c(1:15,j+15,j+16)])
  j=j+2
}

###Divide between the variants present only in the donor and the variants only present in  the recipient
for (i in 1:28){
  assign(paste("df_joint_pair",(i),"_presentDonor",sep=""),
         get(paste("df_joint_pair",(i),sep=""))[which((get(paste("df_joint_pair",(i),sep=""))[,17]=="0" & get(paste("df_joint_pair",(i),sep=""))[,16]!="0") | 
                                                        (get(paste("df_joint_pair",(i),sep=""))[,17]=="1" & (get(paste("df_joint_pair",(i),sep=""))[,17]=="2"))),])
  assign(paste("df_joint_pair",(i),"_presentRecipient",sep=""),
         get(paste("df_joint_pair",(i),sep=""))[which((get(paste("df_joint_pair",(i),sep=""))[,16]=="0" & get(paste("df_joint_pair",(i),sep=""))[,17]!="0") | 
                                                        (get(paste("df_joint_pair",(i),sep=""))[,16]=="1" & (get(paste("df_joint_pair",(i),sep=""))[,17]=="2"))),])
  
}


##Count number of variants per pair that are DIFFERENT
num_presentDonor <- NULL
num_presentRecipient <- NULL
for (i in 1:28) {
  num_presentDonor[i]<- dim(get(paste("df_joint_pair",(i),"_presentDonor",sep="")))[1]
  num_presentRecipient[i]<-dim(get(paste("df_joint_pair",(i),"_presentRecipient",sep="")))[1]
}

id<-seq(16,ncol(df_joint_qc),2)
xx<-cbind(num_presentDonor,num_presentRecipient)
rownames(xx)<-colnames(df_joint_qc)[id]
write.table(xx,file="/Users/Pinedasans/Data/Catalyst/ExomeSeq/Total_variants.txt",row.names = T,col.names = T,sep="\t")

##Save in an Rdata to work further
save(list = ls(all.names = TRUE), file = "/Users/Pinedasans/Data/Catalyst/ExomeSeq/DiffByPairs.RData")

