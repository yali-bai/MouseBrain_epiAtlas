#########    All "our" in the following code refers to Joint Cabernet.
library(Seurat)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""


#Extract our RNA data after integration with zeng
merged.seuratobj.sct=readRDS(paste0(indir,"/integrated_selected_brain_region_of_zeng.top_1000markers.rds"))
our.df <- as.matrix(merged.seuratobj.sct@assays$RNA$counts)[,which(merged.seuratobj.sct$source == "our")]
our.seuratobj <- CreateSeuratObject(our.df)
our.seuratobj  <- NormalizeData(our.seuratobj , verbose = FALSE)
metainfo<-read.csv("../../04data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv",row.names=1)
metainfo<-metainfo[metainfo$RNA_QC==1,]
metainfo<-metainfo[gsub(".*@@_","",colnames(our.seuratobj)),]
# identical(colnames(our.seuratobj),rownames(metainfo))
our.seuratobj$age<-metainfo$old_young
our.seuratobj$lt_twice_class<-metainfo$class_label
our.seuratobj$lt_twice_subclass<-metainfo$subclass_label
saveRDS(our.seuratobj,paste0(outdir,"/Joint_Cabernet_seurat.rds"))
#preprocessing
our.seuratobj <- SCTransform(object = our.seuratobj)
our.seuratobj <- FindVariableFeatures(our.seuratobj , selection.method = "vst", nfeatures = 3000)
our.seuratobj <- RunPCA(our.seuratobj , features = VariableFeatures(object = our.seuratobj ))
our.seuratobj <- RunUMAP(our.seuratobj , dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40)
our.seuratobj <- RunTSNE(our.seuratobj , dims=1:30, dim.embed=3, perplexity=25)
our.seuratobj <- FindNeighbors(our.seuratobj , dims=1:30, reduction="pca")
our.seuratobj <- FindClusters(our.seuratobj)
saveRDS(our.seuratobj,paste0(outdir,"/Joint_Cabernet_seruat_with_cluster_corrected.rds"))



#Filter out low-quality genes
# The filter condition of genes is that the count value is greater than 500 and the number of cells with cpm value greater than 0.5 is greater than 2.
RNA_count<-readRDS(paste0(indir,"/RNA_matrix.rds"))
rowsum<-rowSums(RNA_count)
RNA_filter_count_genes<-rownames(RNA_count)[rowsum>500]

colsum<-colSums(RNA_count)
RNA_cpm<-(RNA_count+0.0001)/colsum*1e6
RNA_filter_cpm_genes<-rownames(RNA_cpm)[rowSums(RNA_cpm>0.5)>=2]

filtered_genes<-intersect(RNA_filter_count_genes,RNA_filter_cpm_genes)
#The crossover gene is the result of RNA_filter_count_genes, which filters out the remaining genes.
write.csv(filtered_genes,"../../output/03.Aging_Mouse/07-aged_DEG/count_cpm_filtered_gene.csv",row.names=F)