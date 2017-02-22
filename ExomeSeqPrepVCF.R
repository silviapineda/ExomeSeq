rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Esome Sequencing Donor/Recipient pairs
###
### CITATION: 
###
### PROCESS: Reading and Treating the vcf files to convert to a data frame with genotypes
###          Reading the annotation file
###           Put all the individuals together in the same data frame
###           
### DESCRIP: ExomeSeq Data from Donors/Recipients pairs 
###           27 pairs + 1 pair with two donors
###         
###
### Author: Silvia Pineda
### Date: June 1, 2016
############################################################################################



#########################
### Installing package ###
#########################
#source("http://bioconductor.org/biocLite.R") 
# biocLite("VariantAnnotation") #install the package
# biocLite('snpStats')
# biocLite("SNPRelate")

#Work Directory
setwd("~/Data/Catalyst/ExomeSeq/VCF/")

library("VariantAnnotation") #load the package
joint.vcf <- readVcf("joint.output.snp.indel.recal.filtered.biallelic.vcf","hg19")

                

library(SNPRelate)
###Convert into a data frame
snpgdsVCF2GDS("joint.output.snp.indel.recal.filtered.biallelic.vcf", "joint.output.snp.indel.recal.filtered.biallelic.gds", method="copy.num.of.ref")
genofile <- snpgdsOpen("joint.output.snp.indel.recal.filtered.biallelic.gds")
for (i in 1:55){
  assign(paste("sample_id_",i,sep=""),read.gdsn(index.gdsn(genofile, "genotype"))[i,])
}

df_joint <- data.frame(
  snp_id = rownames(joint.vcf),
  snp_position = read.gdsn(index.gdsn(genofile, "snp.position")),
  snp_chromosome = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
  snp_allele = read.gdsn(index.gdsn(genofile, "snp.allele"))
)

df_joint <- cbind(df_joint,sample_id_1,sample_id_2,sample_id_3,sample_id_4,sample_id_5,sample_id_6,sample_id_7,sample_id_8,sample_id_9,sample_id_10,
      sample_id_11,sample_id_12,sample_id_13,sample_id_14,sample_id_15,sample_id_16,sample_id_17,sample_id_18,sample_id_19,sample_id_20,
      sample_id_21,sample_id_22,sample_id_23,sample_id_24,sample_id_25,sample_id_26,sample_id_27,sample_id_28,sample_id_29,sample_id_30,
      sample_id_31,sample_id_32,sample_id_33,sample_id_34,sample_id_35,sample_id_36,sample_id_37,sample_id_38,sample_id_39,sample_id_40,
      sample_id_41,sample_id_42,sample_id_43,sample_id_44,sample_id_45,sample_id_46,sample_id_47,sample_id_48,sample_id_49,sample_id_50,
      sample_id_51,sample_id_52,sample_id_53,sample_id_54,sample_id_55)

splitpop <- strsplit(as.character(read.gdsn(index.gdsn(genofile, "sample.id"))),"_")
colnames(df_joint)[5:59] <- unlist(lapply(splitpop, "[", 1))

##Change 0 by 2 and 3 by NA
genotypes <- sapply(5:dim(df_joint)[2],function(x) replace(df_joint[,x],df_joint[,x]=="3","NA"))
genotypes.2 <- sapply(1:dim(genotypes)[2],function(x) replace(genotypes[,x],genotypes[,x]=="2","3"))
genotypes.3 <- sapply(1:dim(genotypes.2)[2],function(x) replace(genotypes.2[,x],genotypes.2[,x]=="0","2"))
genotypes.4 <- sapply(1:dim(genotypes.3)[2],function(x) replace(genotypes.3[,x],genotypes.3[,x]==3,"0"))

df_joint[,5:dim(df_joint)[2]] <- genotypes.4

#Filter by chromosome
df_joint <- df_joint[which(df_joint$snp_chromosome!= "X" & df_joint$snp_chromosome!= "Y" & df_joint$snp_chromosome!= "MT" & df_joint$snp_chromosome!= "GL000191.1" & df_joint$snp_chromosome!= "GL000192.1"  & df_joint$snp_chromosome!= "GL000193.1" 
                     & df_joint$snp_chromosome!= "GL000194.1" & df_joint$snp_chromosome!= "GL000195.1" 
                     & df_joint$snp_chromosome!= "GL000196.1" & df_joint$snp_chromosome!= "GL000197.1" & df_joint$snp_chromosome!= "GL000198.1" & df_joint$snp_chromosome!= "GL000199.1" & df_joint$snp_chromosome!= "GL000201.1" & df_joint$snp_chromosome!= "GL000202.1" 
                     & df_joint$snp_chromosome!= "GL000203.1" & df_joint$snp_chromosome!= "GL000204.1" & df_joint$snp_chromosome!= "GL000205.1" & df_joint$snp_chromosome!= "GL000206.1" & df_joint$snp_chromosome!= "GL000207.1" & df_joint$snp_chromosome!= "GL000208.1" 
                     & df_joint$snp_chromosome!= "GL000209.1"  & df_joint$snp_chromosome!= "GL000210.1" & df_joint$snp_chromosome!= "GL000211.1" & df_joint$snp_chromosome!= "GL000212.1" & df_joint$snp_chromosome!= "GL000213.1" & df_joint$snp_chromosome!= "GL000214.1" 
                     & df_joint$snp_chromosome!="GL000215.1" & df_joint$snp_chromosome!= "GL000216.1" & df_joint$snp_chromosome!= "GL000217.1" & df_joint$snp_chromosome!= "GL000218.1" & df_joint$snp_chromosome!= "GL000219.1" & df_joint$snp_chromosome!= "GL000220.1" 
                     & df_joint$snp_chromosome!="GL000221.1" & df_joint$snp_chromosome!= "GL000222.1" & df_joint$snp_chromosome!= "GL000223.1" & df_joint$snp_chromosome!= "GL000224.1" & df_joint$snp_chromosome!= "GL000225.1" & df_joint$snp_chromosome!= "GL000226.1" 
                     & df_joint$snp_chromosome!="GL000227.1" & df_joint$snp_chromosome!="GL000228.1" & df_joint$snp_chromosome!= "GL000229.1" & df_joint$snp_chromosome!= "GL000230.1" & df_joint$snp_chromosome!= "GL000231.1" & df_joint$snp_chromosome!= "GL000232.1" 
                     & df_joint$snp_chromosome!= "GL000233.1" & df_joint$snp_chromosome!="GL000234.1" & df_joint$snp_chromosome!= "GL000235.1" & df_joint$snp_chromosome!= "GL000236.1" & df_joint$snp_chromosome!= "GL000237.1" & df_joint$snp_chromosome!= "GL000238.1" 
                     & df_joint$snp_chromosome!= "GL000239.1" & df_joint$snp_chromosome!= "GL000240.1" & df_joint$snp_chromosome!="GL000241.1" & df_joint$snp_chromosome!= "GL000242.1" & df_joint$snp_chromosome!= "GL000243.1" & df_joint$snp_chromosome!= "GL000244.1" 
                     & df_joint$snp_chromosome!= "GL000245.1" & df_joint$snp_chromosome!= "GL000246.1" & df_joint$snp_chromosome!= "GL000247.1" & df_joint$snp_chromosome!= "GL000248.1"),]

#n=568,706

###Filther those that are not annotated 
my_annovar<-read.csv("~/Data/Catalyst/ExomeSeq/ANNOVAR/myanno.joint.output.snp.indel.recal.avinput.hg19_multianno.csv")
colnames(df_joint)[2:3]<-c("Start","Chr")
merge_data <- merge(my_annovar,df_joint, by=c("Chr","Start")) #536,061
merge_data_annotated <- merge_data[,c(1:7,9,14:18,44:98,100,99,100)]
colnames(merge_data_annotated)[c(69,71)]<-c("Discovery-270862-R1","Discovery-270862-R2")
df_joint_annotated<-merge_data_annotated

##QC of 95% of non-missingness
exome_variants<-data.matrix(df_joint_annotated[,16:71])
rownames(exome_variants)<-df_joint_annotated$snp_id
exome_variants_qc<-exome_variants[rowSums(is.na(exome_variants) == FALSE) >= 52,] #507,218 that has at least 95% of non-missingness
id.qc<-match(rownames(exome_variants_qc),df_joint_annotated$snp_id)
df_joint_qc<-df_joint_annotated[na.omit(id.qc),]


##Save in the Rdata
save(df_joint_qc,exome_variants_qc,file="ExomeSeqVCF.Rdata")

##To obatin the same list of variants with the 1000G
write.table(df_joint_qc[,1:2],"variantList_vcf.txt",row.names = FALSE,sep="\t")
############################################################################################################################################################

####################
#### Only SNPs #####
####################

###Convert into a data frame
snpgdsVCF2GDS("joint.output.snp.indel.recal.filtered.biallelic.SNPs.vcf", "joint.output.snp.indel.recal.filtered.biallelic.SNPs.gds", method="copy.num.of.ref")
genofile <- snpgdsOpen("joint.output.snp.indel.recal.filtered.biallelic.SNPs.gds")
for (i in 1:55){
  assign(paste("sample_id_",i,sep=""),read.gdsn(index.gdsn(genofile, "genotype"))[i,])
}

df_joint <- data.frame(
  snp_id = rownames(joint.vcf),
  snp_position = read.gdsn(index.gdsn(genofile, "snp.position")),
  snp_chromosome = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
  snp_allele = read.gdsn(index.gdsn(genofile, "snp.allele"))
)

df_joint <- cbind(df_joint,sample_id_1,sample_id_2,sample_id_3,sample_id_4,sample_id_5,sample_id_6,sample_id_7,sample_id_8,sample_id_9,sample_id_10,
                  sample_id_11,sample_id_12,sample_id_13,sample_id_14,sample_id_15,sample_id_16,sample_id_17,sample_id_18,sample_id_19,sample_id_20,
                  sample_id_21,sample_id_22,sample_id_23,sample_id_24,sample_id_25,sample_id_26,sample_id_27,sample_id_28,sample_id_29,sample_id_30,
                  sample_id_31,sample_id_32,sample_id_33,sample_id_34,sample_id_35,sample_id_36,sample_id_37,sample_id_38,sample_id_39,sample_id_40,
                  sample_id_41,sample_id_42,sample_id_43,sample_id_44,sample_id_45,sample_id_46,sample_id_47,sample_id_48,sample_id_49,sample_id_50,
                  sample_id_51,sample_id_52,sample_id_53,sample_id_54,sample_id_55)

splitpop <- strsplit(as.character(read.gdsn(index.gdsn(genofile, "sample.id"))),"_")
colnames(df_joint)[5:59] <- unlist(lapply(splitpop, "[", 1))

##Change 0 by 2 and 3 by NA
genotypes <- sapply(5:dim(df_joint)[2],function(x) replace(df_joint[,x],df_joint[,x]=="3","NA"))
genotypes.2 <- sapply(1:dim(genotypes)[2],function(x) replace(genotypes[,x],genotypes[,x]=="2","3"))
genotypes.3 <- sapply(1:dim(genotypes.2)[2],function(x) replace(genotypes.2[,x],genotypes.2[,x]=="0","2"))
genotypes.4 <- sapply(1:dim(genotypes.3)[2],function(x) replace(genotypes.3[,x],genotypes.3[,x]==3,"0"))

df_joint[,5:dim(df_joint)[2]] <- genotypes.4

#Filter by chromosome
df_joint <- df_joint[which(df_joint$snp_chromosome!= "X" & df_joint$snp_chromosome!= "Y" & df_joint$snp_chromosome!= "MT" & df_joint$snp_chromosome!= "GL000191.1" & df_joint$snp_chromosome!= "GL000192.1"  & df_joint$snp_chromosome!= "GL000193.1" 
                           & df_joint$snp_chromosome!= "GL000194.1" & df_joint$snp_chromosome!= "GL000195.1" 
                           & df_joint$snp_chromosome!= "GL000196.1" & df_joint$snp_chromosome!= "GL000197.1" & df_joint$snp_chromosome!= "GL000198.1" & df_joint$snp_chromosome!= "GL000199.1" & df_joint$snp_chromosome!= "GL000201.1" & df_joint$snp_chromosome!= "GL000202.1" 
                           & df_joint$snp_chromosome!= "GL000203.1" & df_joint$snp_chromosome!= "GL000204.1" & df_joint$snp_chromosome!= "GL000205.1" & df_joint$snp_chromosome!= "GL000206.1" & df_joint$snp_chromosome!= "GL000207.1" & df_joint$snp_chromosome!= "GL000208.1" 
                           & df_joint$snp_chromosome!= "GL000209.1"  & df_joint$snp_chromosome!= "GL000210.1" & df_joint$snp_chromosome!= "GL000211.1" & df_joint$snp_chromosome!= "GL000212.1" & df_joint$snp_chromosome!= "GL000213.1" & df_joint$snp_chromosome!= "GL000214.1" 
                           & df_joint$snp_chromosome!="GL000215.1" & df_joint$snp_chromosome!= "GL000216.1" & df_joint$snp_chromosome!= "GL000217.1" & df_joint$snp_chromosome!= "GL000218.1" & df_joint$snp_chromosome!= "GL000219.1" & df_joint$snp_chromosome!= "GL000220.1" 
                           & df_joint$snp_chromosome!="GL000221.1" & df_joint$snp_chromosome!= "GL000222.1" & df_joint$snp_chromosome!= "GL000223.1" & df_joint$snp_chromosome!= "GL000224.1" & df_joint$snp_chromosome!= "GL000225.1" & df_joint$snp_chromosome!= "GL000226.1" 
                           & df_joint$snp_chromosome!="GL000227.1" & df_joint$snp_chromosome!="GL000228.1" & df_joint$snp_chromosome!= "GL000229.1" & df_joint$snp_chromosome!= "GL000230.1" & df_joint$snp_chromosome!= "GL000231.1" & df_joint$snp_chromosome!= "GL000232.1" 
                           & df_joint$snp_chromosome!= "GL000233.1" & df_joint$snp_chromosome!="GL000234.1" & df_joint$snp_chromosome!= "GL000235.1" & df_joint$snp_chromosome!= "GL000236.1" & df_joint$snp_chromosome!= "GL000237.1" & df_joint$snp_chromosome!= "GL000238.1" 
                           & df_joint$snp_chromosome!= "GL000239.1" & df_joint$snp_chromosome!= "GL000240.1" & df_joint$snp_chromosome!="GL000241.1" & df_joint$snp_chromosome!= "GL000242.1" & df_joint$snp_chromosome!= "GL000243.1" & df_joint$snp_chromosome!= "GL000244.1" 
                           & df_joint$snp_chromosome!= "GL000245.1" & df_joint$snp_chromosome!= "GL000246.1" & df_joint$snp_chromosome!= "GL000247.1" & df_joint$snp_chromosome!= "GL000248.1"),]

#n=

###Filther those that are not annotated 
my_annovar<-read.csv("~/Data/Catalyst/ExomeSeq/ANNOVAR/myanno.joint.output.recal.filtered.biallelic.SNPs.avinput.hg19_multianno.csv")
colnames(df_joint)[2:3]<-c("Start","Chr")
merge_data <- merge(my_annovar,df_joint, by=c("Chr","Start")) #536,061
merge_data_annotated <- merge_data[,c(1:7,9,14:18,44:98,100,99,100)]
colnames(merge_data_annotated)[c(69,71)]<-c("Discovery-270862-R1","Discovery-270862-R2")
df_joint_annotated<-merge_data_annotated

##QC of 95% of non-missingness
exome_variants<-data.matrix(df_joint_annotated[,16:71])
rownames(exome_variants)<-df_joint_annotated$snp_id
exome_variants_qc<-exome_variants[rowSums(is.na(exome_variants) == FALSE) >= 52,] #507,218 that has at least 95% of non-missingness
id.qc<-match(rownames(exome_variants_qc),df_joint_annotated$snp_id)
df_joint_qc<-df_joint_annotated[na.omit(id.qc),]

##Save in the Rdata
save(df_joint_qc,exome_variants_qc,file="ExomeSeqVCF_SNPs.Rdata")



#########################
####Read the 1000G vcf 
#########################
library("VariantAnnotation") #load the package

variants.1000G.vcf <- readVcf("1000G.variants.selected.vcf","hg19")


library(SNPRelate)
###Convert into a data frame
snpgdsVCF2GDS("1000G.variants.selected.vcf", "1000G.variants.selected.vcf.gds", method="copy.num.of.ref")
genofile <- snpgdsOpen("1000G.variants.selected.vcf.gds")
df.1000G <- data.frame(
  snp_id = rownames(variants.1000G.vcf),
  snp_position = read.gdsn(index.gdsn(genofile, "snp.position")),
  snp_chromosome = read.gdsn(index.gdsn(genofile, "snp.chromosome")),
  snp_allele = read.gdsn(index.gdsn(genofile, "snp.allele")),
  snp_annot = read.gdsn(index.gdsn(genofile,  path="snp.annot/filter"))
)


df.1000G <- cbind(df.1000G,t(read.gdsn(index.gdsn(genofile, "genotype"))))
colnames(df.1000G)[6:dim(df.1000G)[2]]<-colnames(variants.1000G.vcf)
save(df.1000G,file="ExomeSeq1000GVCF.Rdata")

write.table(df.1000G,file="df.1000G.txt")
##I have executed the program Mikel has sent to me (silvia-project-assembly-1.0.jar) to chande the 0 for 2 and 2 for 0.
##The name is df.1000G.mod.txt
