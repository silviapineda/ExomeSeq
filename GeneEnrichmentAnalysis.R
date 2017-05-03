rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Annotation
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: August, 2016
############################################################################################

setwd("/Users/Pinedasans/Catalyst/Results/")

#########################
####Read the annotation
#######################
###Read the file with the results
resultsRF<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointRF.txt",header=T,sep="\t")
resultsRF_assocRej<-resultsRF[which(resultsRF$OR>1.4),]

resultFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTestSign.txt",header=T,sep="\t")
resultFisher_assocRej<-resultFisher[which(resultFisher$OR>1),]

riskAMR<-read.table("/Users/Pinedasans/Catalyst/Results/RiskAMR.txt",header=T,sep="\t")
riskCMR<-read.table("/Users/Pinedasans/Catalyst/Results/RiskCMR.txt",header=T,sep="\t")
riskNoRej<-read.table("/Users/Pinedasans/Catalyst/Results/ProtectionNoRej.txt",header=T,sep="\t")

id.overlap<-match(resultFisher_assocRej$Gene.refGene,resultsRF_assocRej$Gene.refGene)
resultFisher_assocRej[which(is.na(id.overlap)==T),]

####Kidney
gene.expr.up.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/results.sign.upregulated.kidney.txt")
gene.expr.each.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/genes.kidney.each.txt")
gene.expr.high.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/gene.kidney.high.expr.txt")
gene.kidney<-union(union(rownames(gene.expr.up.kidney),as.character(gene.expr.each.kidney[,1])),rownames(gene.expr.high.kidney)) #2967

length(na.omit(match(riskAMR[,1],gene.kidney))) #13
length(na.omit(match(riskCMR[,1],gene.kidney))) #3


##Vessels
gene.expr.up.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/results.sign.upregulated.vessels.txt")
gene.expr.each.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/genes.vessels.each.txt")
gene.expr.high.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/gene.vessel.high.expr.txt")
gene.vessels<-union(union(rownames(gene.expr.up.vessels),as.character(gene.expr.each.vessels[,1])),rownames(gene.expr.high.vessels)) #3513

length(na.omit(match(riskAMR[,1],gene.vessels))) #15
length(na.omit(match(riskCMR[,1],gene.vessels))) #6


###Immuno
Dataset1<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/InnateDB_genes1.csv",header=T)
ListGene1<-na.omit(as.character(Dataset1$name))
Dataset2<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/InnateDB_genes2.csv",header=T)
ListGene2<-na.omit(as.character(Dataset2$name))
Dataset3<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/InnateDB_genes3.csv",header=T)
ListGene3<-na.omit(as.character(Dataset3$SourceGeneName))
Dataset4<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/InnateDB_genes4.csv",header=T)
ListGene4<-na.omit(as.character(Dataset4$name))
gene.immuno<-union(union(union(ListGene1,ListGene2),ListGene3),ListGene4) #8745

length(na.omit(match(riskAMR[,1],gene.immuno))) #20
length(na.omit(match(riskCMR[,1],gene.immuno))) #4


##Surface genes
Dataset1<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/peptideatlas.csv",header=T)
ListGene1<-na.omit(as.character(Dataset1$Gene_Symbol))
ListGene1<-ListGene1[which(ListGene1!="")]
Dataset2<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/proteomes.csv",header=T)
ListGene2<-na.omit(as.character(Dataset2$gene_symbol))
ListGene2<-ListGene2[which(ListGene2!="-")]
Dataset3<-read.csv("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/surfaceome.csv",header=T)
ListGene3<-na.omit(as.character(Dataset3$SourceGeneName))

# Dataset12<-intersect(unique(ListGene1), unique(ListGene2))
# Dataset13<-intersect(unique(ListGene1), unique(ListGene3))
# Dataset23<-intersect(unique(ListGene2), unique(ListGene3))
# gene.surface<-unique(c(Dataset12,Dataset13,Dataset23)) #1491

gene.surface<-union(union(ListGene1,ListGene2),ListGene3)  #7341

length(na.omit(match(riskAMR[,1],gene.surface))) #48
length(na.omit(match(riskCMR[,1],gene.surface))) #11


##Enrichment analysis #19725 total coding genes
##The number of genes unique AMR is 72

counts = (matrix(data = c(13, 59, 2786, 16939), nrow = 2)) #p-value = 0.5
chisq.test(counts)

counts = (matrix(data = c(15, 57, 3291, 16434), nrow = 2)) #p-value = 0.5
chisq.test(counts)

counts = (matrix(data = c(20, 52, 8745, 10980), nrow = 2)) #p-value = 0.0003
chisq.test(counts)

counts = (matrix(data = c(48,24, 7341, 12384), nrow = 2)) #p-value = 1.7 *10-6
chisq.test(counts)

##The number of genes unique CMR is 22

counts = (matrix(data = c(3, 19, 2786, 16939), nrow = 2)) #p-value = 0.5
chisq.test(counts)

counts = (matrix(data = c(6, 16, 3291, 16434), nrow = 2)) #p-value = 0.5
chisq.test(counts)

counts = (matrix(data = c(4, 18, 8745, 10980), nrow = 2)) #p-value = 0.0003
chisq.test(counts)

counts = (matrix(data = c(11,11, 7341, 12384), nrow = 2)) #p-value = 1.7 *10-6
chisq.test(counts)

genes_GTEx<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/gencode.v19.annotation_coding_genes.txt",sep=";")
genes_GTEx<-genes_GTEx[,c(1,2,3,4,6)]
colnames(genes_GTEx)<-c("Chr","Start","end","transcipt","gene")

eQTL_sign_artery_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Aorta_eQTLs_sign.txt")
colnames(eQTL_sign_artery_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_artery_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_artery_aorta2<-cbind(eQTL_sign_artery_aorta,genes_GTEx[id.transript,])
eQTL_sign_artery_aorta_coding<-eQTL_sign_artery_aorta2[which(is.na(id.transript)==F),]
colnames(eQTL_sign_artery_aorta_coding)[8:10]<-c("gene_chr","gene_start","gene_end")

eQTL_sign_Coronary_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Coronary_eQTLs_sign.txt")
colnames(eQTL_sign_Coronary_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Coronary_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_Coronary_aorta2<-cbind(eQTL_sign_Coronary_aorta,genes_GTEx[id.transript,])
eQTL_sign_Coronary_aorta_coding<-eQTL_sign_Coronary_aorta2[which(is.na(id.transript)==F),]
colnames(eQTL_sign_Coronary_aorta_coding)[8:10]<-c("gene_chr","gene_start","gene_end")

eQTL_sign_Tibial_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Tibial_eQTLs_sign.txt")
colnames(eQTL_sign_Tibial_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Tibial_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_Tibial_aorta2<-cbind(eQTL_sign_Tibial_aorta,genes_GTEx[id.transript,])
eQTL_sign_Tibial_aorta_coding<-eQTL_sign_Tibial_aorta2[which(is.na(id.transript)==F),]
colnames(eQTL_sign_Tibial_aorta_coding)[8:10]<-c("gene_chr","gene_start","gene_end")

eQTL_bloodVessels<-rbind(eQTL_sign_artery_aorta_coding,eQTL_sign_Coronary_aorta_coding,eQTL_sign_Tibial_aorta_coding)

merge_bloodVessels<-merge(resultFisher_assocRej,eQTL_bloodVessels,by=c("Chr","Start"))
write.table(merge_bloodVessels,"/Users/Pinedasans/Catalyst/Results/eQTL_bloodVessels.txt",sep="\t",row.names = F)

##Delete manually the repeated genes by transcript and leave the most signifcant transcript
merge_bloodVessels<-read.table("/Users/Pinedasans/Catalyst/Results/eQTL_bloodVessels.txt",sep="\t",header = T)

length(na.omit(match(merge_bloodVessels$gene,gene.kidney))) #9
length(na.omit(match(merge_bloodVessels$gene,gene.vessels))) #9
length(na.omit(match(merge_bloodVessels$gene,gene.immuno))) #5
length(na.omit(match(merge_bloodVessels$gene,gene.surface))) #4

##Enrichment analysis #19725 total coding genes
counts = (matrix(data = c(9, 16, 2786, 16939), nrow = 2)) #p-value= 0.004
chisq.test(counts)

counts = (matrix(data = c(9, 16, 3291, 16434), nrow = 2)) #p-value = 0.02
chisq.test(counts)

counts = (matrix(data = c(5, 19, 4677, 15058), nrow = 2)) #p-value = 0.9
chisq.test(counts)

counts = (matrix(data = c(4,20, 3845, 15880), nrow = 2)) #p-value = 0.9
chisq.test(counts)

eQTL_sign_Whole_blood<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Whole_Blood_eQTLs_sign.txt")
colnames(eQTL_sign_Whole_blood)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Whole_blood$transcript,genes_GTEx$transcipt)
eQTL_sign_Whole_blood2<-cbind(eQTL_sign_Whole_blood,genes_GTEx[id.transript,])
eQTL_sign_Whole_blood_coding<-eQTL_sign_Whole_blood2[which(is.na(id.transript)==F),]
colnames(eQTL_sign_Whole_blood_coding)[8:10]<-c("gene_chr","gene_start","gene_end")

merge_wholeBlood<-merge(resultFisher_assocRej,eQTL_sign_Whole_blood_coding,by=c("Chr","Start")) #17variants
write.table(merge_wholeBlood,"/Users/Pinedasans/Catalyst/Results/eQTL_wholeBlood.txt",sep="\t",row.names = F)

merge_bloodVessels<-read.table("/Users/Pinedasans/Catalyst/Results/eQTL_bloodVessels.txt",sep="\t",header = T)

length(na.omit(match(merge_wholeBlood$gene,gene.kidney))) #7
length(na.omit(match(merge_wholeBlood$gene,gene.vessels))) #7
length(na.omit(match(merge_wholeBlood$gene,gene.immuno))) #3
length(na.omit(match(merge_wholeBlood$gene,gene.surface))) #2

##Enrichment analysis #19725 total coding genes
counts = (matrix(data = c(7, 7, 2786, 16939), nrow = 2)) #p-value = 0.0005
chisq.test(counts)

counts = (matrix(data = c(7, 7, 3291, 16434), nrow = 2)) #p-alue = 0.002
chisq.test(counts)

counts = (matrix(data = c(3, 11, 4677, 15058), nrow = 2)) #p-value=0.999
chisq.test(counts)

counts = (matrix(data = c(2,12, 3845, 15880), nrow = 2)) #p-value=0.8
chisq.test(counts)







###Within each list - genes
##List Endpoint
library("VennDiagram")
area1=length(which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES"))
area2=length(which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES"))
area3=length(which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES"))
area4=length(which(annotated.variants.assoc.endpoint.AMR$Surface=="YES"))
n12=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")]))
n13=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]))
n14=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n23=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]))
n24=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n34=length(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")],
                     annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n123=length(intersect(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                                annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")]),
                      annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]))
n124=length(intersect(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                                annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")]),
                      annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n134=length(intersect(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                                annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]),
                      annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n234=length(intersect(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")],
                                annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]),
                      annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))
n1234=length(intersect(intersect(intersect(annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Kidney=="YES")],
                                           annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Vessels=="YES")]),
                                 annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Immuno=="YES")]),
                       annotated.variants.assoc.endpoint.AMR$Gene.refGene[which(annotated.variants.assoc.endpoint.AMR$Surface=="YES")]))

tiff("/Users/Pinedasans/Documents/Catalyst/Results_V3/VennDiagramEndpointAMR.tiff",width=95,height=95,units="mm",res=300,compression=c("lzw"))

draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
               n34, n123, n124, n134, n234, n1234,euler.d=F,scaled=F,
               category = c("Kidney (23)","Blood Vessels (33)","Immune (49)","Surface (25)" ),
               lty = "blank", fill=rainbow(4), alpha = rep(0.5, 4), cex = 0.7,lwd = rep(1, 4),
               cat.cex=rep(0.5,4))

dev.off()








