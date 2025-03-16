#########    All "our" in the following code refers to Joint Cabernet.

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

#############  Sort out the metainfo
# The mean, cv and na ratio of all genes in RNA were calculated.
library(data.table)
library(MuDataSeurat)
library(dplyr)


our_RNA<-readRDS(paste0(indir,"/our_RNA_annotated_latest.rds"))
our_RNA_count<-our_RNA@assays$RNA$data
data<-data.frame(gene_id=rownames(our_RNA_count))
data$RNA.NA.ratio<-apply(our_RNA_count,1,function(x){
  na.ratio<-sum(x==0)/length(x)
  return(na.ratio)
})
data$RNA.mean<-apply(our_RNA_count,1,function(x){
  mean<-mean(x)
  return(mean)
})
data$RNA.cv<-apply(our_RNA_count,1,function(x){
  cv<-sd(x)/mean(x)
  return(cv)
})
write.csv(data,"../../../output/01-youth/02-correlation_calculation/RNA_metainfo_data.csv",row.names = F)

#Generate metainfo of all genes on genebody, including gene id, gene name, chromosome location, gene length, CpG number, etc
gene_bed<-read.table("../input/Gene.output.bed.gz", header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("chromosome", "start", "end","Cpg_number"))
gene_location<-read.csv("../input/gene_location.csv")
gene_location$start<-as.numeric(gene_location$start)-1
colnames(gene_location)[1:2]<-c("gene_id","chromosome")
gene_metainfo<-merge(gene_location,gene_bed)
gene_metainfo$Cpg_number<-gsub("CpG:_","",gene_metainfo$Cpg_number)
gene_metainfo$Gene_length<-gene_metainfo$end-gene_metainfo$start+1
gene_meta_last<-read.csv("../input/gene_genelength_cpgnum.csv")
gene_metainfo<-merge(gene_metainfo,gene_meta_last[,1:2])%>%
mutate(var_dim=rep("genebody",nrow(.)))
gene_metainfo<-gene_metainfo[!duplicated(gene_metainfo),]
write.csv(gene_metainfo,"../../../output/01-youth/02-correlation_calculation/gene_metainfo.csv",row.names=F)



