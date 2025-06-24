# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

#############  Sort out the metainfo
# The mean, cv and na ratio of all genes in RNA were calculated.
library(data.table)
library(MuDataSeurat)
library(dplyr)


Joint_Cabernet_RNA<-readRDS(paste0(indir,"/Joint_Cabernet_RNA_annotated_latest.rds"))
Joint_Cabernet_RNA_count<-Joint_Cabernet_RNA@assays$RNA$data
data<-data.frame(gene_id=rownames(Joint_Cabernet_RNA_count))
data$RNA.NA.ratio<-apply(Joint_Cabernet_RNA_count,1,function(x){
  na.ratio<-sum(x==0)/length(x)
  return(na.ratio)
})
data$RNA.mean<-apply(Joint_Cabernet_RNA_count,1,function(x){
  mean<-mean(x)
  return(mean)
})
data$RNA.cv<-apply(Joint_Cabernet_RNA_count,1,function(x){
  cv<-sd(x)/mean(x)
  return(cv)
})
write.csv(data,"../../../output/01.Young_Mouse/02-correlation_calculation/RNA_metainfo_data.csv",row.names = F)

#Generate metainfo of all genes on genebody, including gene id, gene name, chromosome location, gene length, CpG number, etc
gene_bed_cpg<-read.table("../../../04.data/01.ref/Gene.output.bed.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chr", "start", "end","Cpg_number"))
gene_bed<-read.csv("../../../04.data/01.ref/mm10.genes_duplicated.bed")
gene_bed$start<-as.numeric(gene_bed$start)-1
gene_metainfo<-merge(gene_bed,gene_bed_cpg)
gene_metainfo$Cpg_number<-gsub("CpG:_","",gene_metainfo$Cpg_number)
gene_metainfo$Gene_length<-gene_metainfo$end-gene_metainfo$start+1
gene_meta_last<-read.csv("../../../04.data/gene_genelength_cpgnum.csv")
gene_metainfo<-gene_metainfo%>%
  mutate(var_dim=rep("genebody",nrow(.)))%>%
  .[!duplicated(.),]
write.csv(gene_metainfo,"../../../output/01.Young_Mouse/02-correlation_calculation/gene_metainfo.csv",row.names=F)



