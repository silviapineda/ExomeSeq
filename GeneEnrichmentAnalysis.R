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

resultFisher<-read.table("/Users/Pinedasans/Catalyst/Results/ResultsEndpointFisherTest.txt",header=T,sep="\t")
resultFisher_assocRej<-resultFisher[which(resultFisher$OR>1.4),]

id.overlap<-match(resultFisher_assocRej$Gene.refGene,resultsRF_assocRej$Gene.refGene)
resultFisher_assocRej[which(is.na(id.overlap)==T),]

####Kidney
gene.expr.up.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/results.sign.upregulated.kidney.txt")
gene.expr.each.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/genes.kidney.each.txt")
gene.expr.high.kidney<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/gene.kidney.high.expr.txt")
gene.kidney<-union(union(rownames(gene.expr.up.kidney),as.character(gene.expr.each.kidney[,1])),rownames(gene.expr.high.kidney)) #2967

length(na.omit(match(resultsRF_assocRej$Gene.refGene,gene.kidney))) #11


##Vessels
gene.expr.up.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/results.sign.upregulated.vessels.txt")
gene.expr.each.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/genes.vessels.each.txt")
gene.expr.high.vessels<-read.table("/Users/Pinedasans/Catalyst/Data/GeneEnrichment/gene.vessel.high.expr.txt")
gene.vessels<-union(union(rownames(gene.expr.up.vessels),as.character(gene.expr.each.vessels[,1])),rownames(gene.expr.high.vessels)) #3513

length(na.omit(match(resultsRF_assocRej$Gene.refGene,gene.vessels))) #11


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

length(na.omit(match(resultsRF_assocRej$Gene.refGene,gene.immuno))) #18


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

length(na.omit(match(resultsRF_assocRej$Gene.refGene,gene.surface))) #3 #22


##Enrichment analysis #19725 total coding genes
counts = (matrix(data = c(13, 57, 2786, 16939), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(15, 55, 3291, 16434), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(17, 53, 8745, 10980), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(48,22, 7341, 12384), nrow = 2))
chisq.test(counts)


genes_GTEx<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/gencode.v19.annotation_coding_genes.txt",sep=";")
genes_GTEx<-genes_GTEx[,c(1,2,3,4,6)]
colnames(genes_GTEx)<-c("Chr","Start","end","transcipt","gene")

eQTL_sign_artery_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Aorta_eQTLs_sign.txt")
colnames(eQTL_sign_artery_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_artery_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_artery_aorta2<-cbind(eQTL_sign_artery_aorta,genes_GTEx$gene[id.transript])
eQTL_sign_artery_aorta_coding<-eQTL_sign_artery_aorta2[which(is.na(eQTL_sign_artery_aorta2$`genes_GTEx$gene[id.transript]`)==F),]
colnames(eQTL_sign_artery_aorta_coding)[8]<-"gene"

eQTL_sign_Coronary_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Coronary_eQTLs_sign.txt")
colnames(eQTL_sign_Coronary_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Coronary_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_Coronary_aorta2<-cbind(eQTL_sign_Coronary_aorta,genes_GTEx$gene[id.transript])
eQTL_sign_Coronary_aorta_coding<-eQTL_sign_Coronary_aorta2[which(is.na(eQTL_sign_Coronary_aorta2$`genes_GTEx$gene[id.transript]`)==F),]
colnames(eQTL_sign_Coronary_aorta_coding)[8]<-"gene"

eQTL_sign_Tibial_aorta<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Artery_Tibial_eQTLs_sign.txt")
colnames(eQTL_sign_Tibial_aorta)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Tibial_aorta$transcript,genes_GTEx$transcipt)
eQTL_sign_Tibial_aorta2<-cbind(eQTL_sign_Tibial_aorta,genes_GTEx$gene[id.transript])
eQTL_sign_Tibial_aorta_coding<-eQTL_sign_Tibial_aorta2[which(is.na(eQTL_sign_Tibial_aorta2$`genes_GTEx$gene[id.transript]`)==F),]
colnames(eQTL_sign_Tibial_aorta_coding)[8]<-"gene"

eQTL_bloodVessels<-rbind(eQTL_sign_artery_aorta_coding,eQTL_sign_Coronary_aorta_coding,eQTL_sign_Tibial_aorta_coding)

merge_bloodVessels<-merge(resultsRF_assocRej,eQTL_bloodVessels,by=c("Chr","Start"))

length(na.omit(match(merge_bloodVessels$Gene.refGene,gene.kidney))) #24
length(na.omit(match(merge_bloodVessels$Gene.refGene,gene.vessels))) #22
length(na.omit(match(merge_bloodVessels$Gene.refGene,gene.immuno))) #9
length(na.omit(match(merge_bloodVessels$Gene.refGene,gene.surface))) #6

##Enrichment analysis #19725 total coding genes
counts = (matrix(data = c(24, 24, 2786, 16939), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(22, 26, 3291, 16434), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(9, 39, 4677, 15058), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(6,42, 3845, 15880), nrow = 2))
chisq.test(counts)

eQTL_sign_Whole_blood<-read.table("/Users/Pinedasans/Catalyst/Data/GTEX/Whole_Blood_eQTLs_sign.txt")
colnames(eQTL_sign_Whole_blood)<-c("Chr","Start","transcript","V4","p-value","Beta","p-vale-adj")
id.transript<-match(eQTL_sign_Whole_blood$transcript,genes_GTEx$transcipt)
eQTL_sign_Whole_blood2<-cbind(eQTL_sign_Whole_blood,genes_GTEx$gene[id.transript])
eQTL_sign_Whole_blood_coding<-eQTL_sign_Whole_blood2[which(is.na(eQTL_sign_Whole_blood2$`genes_GTEx$gene[id.transript]`)==F),]
colnames(eQTL_sign_Whole_blood_coding)[8]<-"gene"

merge_wholeBlood<-merge(resultsRF_assocRej,eQTL_sign_Whole_blood_coding,by=c("Chr","Start")) #17variants

length(na.omit(match(merge_wholeBlood$Gene.refGene,gene.kidney))) #15
length(na.omit(match(merge_wholeBlood$Gene.refGene,gene.vessels))) #12
length(na.omit(match(merge_wholeBlood$Gene.refGene,gene.immuno))) #2
length(na.omit(match(merge_wholeBlood$Gene.refGene,gene.surface))) #5

##Enrichment analysis #19725 total coding genes
counts = (matrix(data = c(15, 10, 2786, 16939), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(12, 13, 3291, 16434), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(2, 25, 4677, 15058), nrow = 2))
chisq.test(counts)

counts = (matrix(data = c(5,25, 3845, 15880), nrow = 2))
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








