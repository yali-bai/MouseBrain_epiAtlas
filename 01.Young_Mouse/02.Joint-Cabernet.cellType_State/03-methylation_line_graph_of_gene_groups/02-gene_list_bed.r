#########    All "our" in the following code refers to Joint Cabernet.
#Generate bed files for gene groups
library(dplyr)
library(Seurat)
library(data.table)
library(MuDataSeurat)
library(tidyverse)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
###########################Obtain criteria for dividing gene groups
our_zeng_union<-readRDS(paste0(indir,"/seuratObj.rds"))
zeng_RNA<-subset(our_zeng_union,subset=group=="zeng")
zeng_RNA_count<-zeng_RNA@assays$RNA$data
zeng_RNA_count<-zeng_RNA_count[rownames(zeng_unique@assays$RNA),]
count_value<- as.vector(zeng_RNA_count)  #变成向量
#Take a non-zero value and find the median
count_value_no_0<-count_value[count_value!=0]
quantiles <- quantile(count_value_no_0, probs = 0.5)


########## Grouping genes and extracting genes from bed files
zeng_RNA_subclass<-read.csv("../../../input/01-youth/zeng_subclass_mean_dat_final.csv",row.names = 1)
bed_data <- read.table("../../../input/reference_genome/gencode.vM18.annotation.gene.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chromosome", "start", "end","strand","gene_id","gene_type","gene_name"))


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

  write.table(group1_bed, file = paste0(outdir,"/",subclass_name,"_group1_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  write.table(group2_bed, file = paste0(outdir,"/",subclass_name,"_group2_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  write.table(group3_bed, file = paste0(outdir,"/",subclass_name,"_group3_bed.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)  
  
}



