library(reshape2)
library(dplyr)
library(Seurat)
library(data.table)
library(MuDataSeurat)
library(tidyverse)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
########################   subclass       ##########################
#RNA
# classes<-c("class","subclass","three_class")
classes<-c("subclass")
age<-c("old")
for(cla in classes){
  for(ag in age){
    print(paste(ag,cla))
    RNA<-read.csv(paste0(indir,"/",cla,"/RNA_",cla,"_mean_dat_final.csv"),row.names = 1)
    RNA$gene<-rownames(RNA)
    RNA_expr<-melt(RNA,id.vars = "gene")
    write.csv(RNA_expr,paste0(outdir,"/",cla,"/RNA_",cla,"_expr.csv"),row.names=F)
  }
}

#DNA
datatype<-c("5mC","5hmC","true_5mC")
var_dim<-c("genebody","promoter")
mctype<-c("CH","CG")
# classes<-c("class","subclass","three_class")
classes<-c("subclass")
age<-c("old","young")
for(cla in classes){
  for(ag in age){
  for(k in var_dim){
    for(j in mctype){
      for (p in datatype) {
        print(paste(k,j,p,ag,cla))
          mc<-read.csv(paste0(indir,"/",cla,"/",p,"_",j,"_",k,"_",cla,"_mean_dat_final.csv"),row.names = 1)
          mc$gene<-rownames(mc)
          mc_expr<-melt(mc,id.vars = "gene")
          write.csv(mc_expr,paste0(outdir,"/",cla,"/",p,"_",j,"_",k,"_",cla,"_expr.csv"),row.names=F)
      }
    }
  }
}
}



