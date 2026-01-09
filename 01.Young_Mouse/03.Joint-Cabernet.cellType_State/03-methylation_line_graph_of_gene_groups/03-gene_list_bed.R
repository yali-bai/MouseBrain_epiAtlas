#Generate bed files for gene groups
library(dplyr)
library(Seurat)
library(data.table)
library(MuDataSeurat)
library(tidyverse)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir="gene_list_bed"
###########################Obtain criteria for dividing gene groups
Joint_Cabernet_zeng_union<-readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/integration_Joint_Cabernet_and_Zeng.with_celltype.rds")
zeng_unique = readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/rds/Zeng_seurat.rds")
zeng_RNA<-subset(Joint_Cabernet_zeng_union,subset=group=="Zeng")
zeng_RNA_count<-zeng_RNA@assays$RNA$data
zeng_RNA_count<-zeng_RNA_count[rownames(zeng_unique@assays$RNA),]
count_value<- as.vector(zeng_RNA_count)  
#Take a non-zero value and find the median
count_value_no_0<-count_value[count_value!=0]
quantiles <- quantile(count_value_no_0, probs = 0.5)


########## Grouping genes and extracting genes from bed files
zeng_RNA_subclass<-read.csv("zeng_subclass_mean_dat_final.csv",row.names = 1)
bed_data <- read.table("../../../03.data/01.ref/mm10.genes.bed", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
outdir = "gene_list_bed"

for(i in 1:ncol(zeng_RNA_subclass)){
  cat(i,'\n')
  data<-zeng_RNA_subclass[,i]%>%as.data.frame()%>%`colnames<-`("mean")%>%mutate(gene_id=rownames(zeng_RNA_subclass))
  subclass_name<-colnames(zeng_RNA_subclass)[i]
  group1<-data[data$mean==0,]
  group2<-data[data$mean>0&data$mean<=quantiles,]
  group3<-data[data$mean>quantiles,]
  
  group1_bed<-bed_data[bed_data$gene_id%in%group1$gene_id,]
  group2_bed<-bed_data[bed_data$gene_id%in%group2$gene_id,]
  group3_bed<-bed_data[bed_data$gene_id%in%group3$gene_id,]
  group1_bed <- group1_bed[,c(1:5,7,6)]
  group2_bed <- group2_bed[,c(1:5,7,6)]
  group3_bed <- group3_bed[,c(1:5,7,6)]

  write.table(group1_bed, file = paste0(outdir,"/",subclass_name,"_group1_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  write.table(group2_bed, file = paste0(outdir,"/",subclass_name,"_group2_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  write.table(group3_bed, file = paste0(outdir,"/",subclass_name,"_group3_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  
}



