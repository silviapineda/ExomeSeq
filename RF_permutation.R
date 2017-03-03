setwd("/home/pinedasans/ExomeSeq/RandomForest/")
load("DataRF_preselected.Rdata")
non.list<- seq(1,56,2) ##Donors

##permutation test
library("VSURF")
for (i in 1:10){
  print(i)    
  perm<-sample(demographics$phenotype[non.list])
  set.seed(2)
  assign(paste("fit",i,sep=""),VSURF(x = dataset, y=perm,parallel = TRUE,ncores=30))
}
save(list = ls(all.names = TRUE),file="ResultsRF_preselected.Rdata")