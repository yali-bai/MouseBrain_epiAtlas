#########    All "our" in the following code refers to Joint Cabernet.
##### 01.import packages #####
library(Seurat)
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library"))
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)
library(getopt)
library(data.table)
library(getopt)
library(Seurat)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.define parameters #####
arguments = matrix(c(
  'help', 'h', 0, "logical", "",
  'num', 'n', 1, "numeric", ""
), byrow=TRUE, ncol=5)
args = getopt(arguments)


if (!is.null(args$help) || is.null(args$num)) {
  cat(paste(getopt(arguments, usage = T), "\n"))
  q()
}

if (!dir.exists(paste0(indir,"/cluster_plot.integrated_by_top_1000_markers"))) {  
    # 如果文件夹不存在，则创建它  
    dir.create(paste0(indir,"/cluster_plot.integrated_by_top_1000_markers"))  
}

##### 03.set working path #####
# setwd("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker")

##### 04.prepare matrix #####
## our matrix old mouse ##
#our_RNA.df <- readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/raw_data/our/old_mouse.raw_count.without_QC_filter.rds")
our_RNA.df <- readRDS("../../input/03-aging/old_mouse.raw_count.without_QC_filter.rds")
#QC.df = fread("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/raw_data/our/old_mouse.RNA.total.stat.txt",header=TRUE,sep="\t",data.table = FALSE)
QC.df = fread("../../input/03-aging/old_mouse.RNA.total.stat.txt",header=TRUE,sep="\t",data.table = FALSE)
QC_subset.df= subset(QC.df,QC == 1) # QC for old mice
our_RNA_subset = our_RNA.df[,match(QC_subset.df$SampleID,colnames(our_RNA.df))]
print(dim(our_RNA_subset)[2])

## Zeng 10X RNA samples and Joint-Cabernet young mice samples ##
merged_samples.df = readRDS(paste0(indir,"/Integrated_RNA_annotated_latest.rds"))

## extract common gene ##
common_gene <-intersect(rownames(our_RNA.df),rownames(merged_samples.df@assays$RNA$counts))
length(common_gene)

## combine matrix ##
zeng.df <- as.matrix(merged_samples.df@assays$RNA$counts[common_gene,])
our_RNA.df <- our_RNA_subset[common_gene,]

print(paste0("zeng data final cell number is ",length(which(merged_samples.df@meta.data$group == "zeng")),"; and our data final cell number is ",length(which(merged_samples.df@meta.data$group == "our"))+dim(our_RNA_subset)[2]))
#which(! rownames(zeng.df) == rownames(our_RNA.df))
merged.df <- cbind(zeng.df,our_RNA.df)
rownames(merged.df) = unlist(lapply(rownames(merged.df), function(x) strsplit(x,"\\.")[[1]][1]))

merged.seuratobj <- CreateSeuratObject(merged.df)
merged.seuratobj$source <- c(merged_samples.df@meta.data$group,rep("our",dim(our_RNA.df)[2])) # define source
rm(merged.df)

# merged.seuratobj = readRDS("merged.seuratobj.rds")
##### 05.integration #####
## step 1. define integration features which is top1000 subclass markers of integration of Zeng 10X RNA samples and Joint-Cabernet young mice samples ##
#marker_gene = readRDS("Integrated_our_zeng_markerGenes.rds") 
marker_gene = readRDS("../../input/03-aging/Integrated_our_zeng_markerGenes.rds") # top1000 subclass markers data
head(marker_gene)
marker_gene$uniq = unlist(lapply(marker_gene$gene, function(x) strsplit(x,"\\.")[[1]][1]))
marker_gene %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        dplyr::filter(p_val_adj < 0.05) %>%
        slice_head(n = args$num) %>%
        ungroup() -> top
print(paste0("the number of gene is ",length(unique(top$uniq))))

merged.seuratobj.list <- SplitObject(merged.seuratobj, split.by = "source")
merged.seuratobj.list <- lapply(X = merged.seuratobj.list, FUN = SCTransform, method = "glmGamPoi")
testfeatures <- SelectIntegrationFeatures(object.list = merged.seuratobj.list)
merged.seuratobj.list <- PrepSCTIntegration(object.list = merged.seuratobj.list, anchor.features = unique(intersect(testfeatures, top$uniq))) 
merged.seuratobj.list <- lapply(X = merged.seuratobj.list, FUN = RunPCA, features = unique(intersect(testfeatures, top$uniq)))
merged.seuratobj.anchors <- FindIntegrationAnchors(object.list = merged.seuratobj.list, normalization.method = "SCT",
                                         anchor.features = unique(intersect(testfeatures, top$uniq)), dims = 1:30, reduction = "rpca", k.anchor = 20)
merged.seuratobj.sct <- IntegrateData(anchorset = merged.seuratobj.anchors, normalization.method = "SCT", dims = 1:30)

## step 2. Dimensionality reduction and clustering##
DefaultAssay(merged.seuratobj.sct) <- "integrated"
merged.seuratobj.sct <- RunPCA(merged.seuratobj.sct, verbose = FALSE)
merged.seuratobj.sct <- RunUMAP(merged.seuratobj.sct, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40)
merged.seuratobj.sct <- RunTSNE(merged.seuratobj.sct, dims=1:30, dim.embed=3, perplexity=25)
merged.seuratobj.sct <- FindNeighbors(merged.seuratobj.sct, dims=1:30, reduction="pca")
merged.seuratobj.sct <- FindClusters(merged.seuratobj.sct)

## step 3. normalize##
DefaultAssay(merged.seuratobj.sct) <- "RNA"
merged.seuratobj.sct <- NormalizeData(merged.seuratobj.sct, verbose = FALSE)


## step 4. integration ##
merged.seuratobj.sct <- JoinLayers(merged.seuratobj.sct)

## save integration result ##
#saveRDS(merged.seuratobj.sct,file=paste0("integrated_selected_brain_region_of_zeng.top_1000markers.rds"))


##### 06. label transfer once with Zeng 10X #####
## read Zeng 10X metainfo in ##
#meta = readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/selected_region/01-zeng_v3_metadata_downsample_1000.rds")
meta = readRDS("../../input/03-aging/01-zeng_v3_metadata_downsample_1000.rds")
meta = as.data.frame(meta)
rownames(meta) = meta$cell_new

## class label transfer ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$seurat_clusters
levels(merged.seuratobj.sct)
result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(merged.seuratobj.sct)) {
  cell_indices <- WhichCells(merged.seuratobj.sct, idents = ident_value)
  table_values <- table(meta[merged.seuratobj.sct@meta.data[cell_indices,"sample"], "class_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
}
print(result_df)

merged.seuratobj.sct$class_label = "NA"
for(i in 1:nrow(result_df)){
  merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$seurat_clusters == result_df$idents[i]),'class_label'] <- result_df$max_variable[i]
}
head(merged.seuratobj.sct@meta.data)


## subclass label transfer ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$seurat_clusters
result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(merged.seuratobj.sct)) {
  cell_indices <- WhichCells(merged.seuratobj.sct, idents = ident_value)
  table_values <- table(meta[merged.seuratobj.sct@meta.data[cell_indices,"sample"], "subclass_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
}
print(result_df)

merged.seuratobj.sct$subclass_label = "NA"
for(i in 1:nrow(result_df)){
  merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$seurat_clusters == result_df$idents[i]),'subclass_label'] <- result_df$max_variable[i]
}
head(merged.seuratobj.sct@meta.data$subclass_label)

##### 06. label transfer twice by Joint-Cabernet cells only #####
our.df <- as.matrix(merged.seuratobj.sct@assays$RNA$counts)[,which(merged.seuratobj.sct$source == "our")]
our.seuratobj <- CreateSeuratObject(our.df)
our.seuratobj <- SCTransform(object = our.seuratobj )
our.seuratobj <- FindVariableFeatures(our.seuratobj , selection.method = "vst", nfeatures = 3000)
our.seuratobj <- RunPCA(our.seuratobj , features = VariableFeatures(object = our.seuratobj ))
our.seuratobj <- RunUMAP(our.seuratobj , dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40)
our.seuratobj <- RunTSNE(our.seuratobj , dims=1:30, dim.embed=3, perplexity=25)
our.seuratobj <- FindNeighbors(our.seuratobj , dims=1:30, reduction="pca")
our.seuratobj <- FindClusters(our.seuratobj) 

## class label transfer twice ##
our.seuratobj$lt_once_class <- merged.seuratobj.sct@meta.data[rownames(our.seuratobj@meta.data),]$class_label # assign class label of integration with Zeng 10X
result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
Idents(our.seuratobj) <- our.seuratobj@meta.data$seurat_clusters
for (ident_value in levels(our.seuratobj)) {
  cell_indices <- WhichCells(our.seuratobj, idents = ident_value)
  table_values <- table(our.seuratobj@meta.data[cell_indices, "lt_once_class"])
  max_variable <- names(table_values)[which.max(table_values)]
  result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
}
print(result_df)

our.seuratobj@meta.data$lt_twice_class = "NA"
for(i in 1:nrow(result_df)){
  our.seuratobj@meta.data[which(our.seuratobj@meta.data$seurat_clusters == result_df$idents[i]),'lt_twice_class'] <- result_df$max_variable[i]
}
head(our.seuratobj@meta.data)

## subclass label transfer twice ##
our.seuratobj$lt_once_subclass <- merged.seuratobj.sct@meta.data[rownames(our.seuratobj@meta.data),]$subclass_label # assign subclass label of integration with Zeng 10X
result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
Idents(our.seuratobj) <- our.seuratobj@meta.data$seurat_clusters
for (ident_value in levels(our.seuratobj)) {
  cell_indices <- WhichCells(our.seuratobj, idents = ident_value)
  table_values <- table(our.seuratobj@meta.data[cell_indices, "lt_once_subclass"])
  max_variable <- names(table_values)[which.max(table_values)]
  result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
}
print(result_df)

our.seuratobj@meta.data$lt_twice_subclass = "NA"
for(i in 1:nrow(result_df)){
  our.seuratobj@meta.data[which(our.seuratobj@meta.data$seurat_clusters == result_df$idents[i]),'lt_twice_subclass'] <- result_df$max_variable[i]
}
head(our.seuratobj@meta.data)

##### 07. add other information #####
## age ##
merged.seuratobj.sct$age = NA
QC.df = fread("../../input/03-aging/old_mouse.RNA.total.stat.txt",header=TRUE,sep="\t",data.table = FALSE)
merged.seuratobj.sct$age[match(QC.df$SampleID,rownames(merged.seuratobj.sct@meta.data))] = "old"
merged.seuratobj.sct$age[which(is.na(merged.seuratobj.sct$age))] = "young"
merged.seuratobj.sct$age[which(merged.seuratobj.sct$source == "zeng")] = NA

## class, subclass label of label transfer twice ##
merged.seuratobj.sct$labeltransfer_twice_class = NA
merged.seuratobj.sct@meta.data[rownames(our.seuratobj@meta.data),"labeltransfer_twice_class"]= our.seuratobj$lt_twice_class
merged.seuratobj.sct$labeltransfer_twice_subclass = NA
merged.seuratobj.sct@meta.data[rownames(our.seuratobj@meta.data),"labeltransfer_twice_subclass"]= our.seuratobj$lt_twice_subclass


##### 08. plot result of integration Joint-Cabernet samples and Zeng 10X #####
#inte.col = readRDS("/share/analysisdata/Methyl/public/TSO/yangfa/analysis/RNA_latest/01-Cluster_analysis/output/01-class_label_our.col.rds")
#Idents(merged.seuratobj.sct) = merged.seuratobj.sct$class_label
#setdiff(unique(merged.seuratobj.sct$class_label),names(inte.col))
#inte.col = c(inte.col,c("CNU-MGE GABA"="#FEE500","Vascular"='#00a6ac','zeng'= 'lightgrey',"our"='lightgrey'))

## major class ##
#inte.col = readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/color.majorclass.rds")
inte.col = readRDS("../../input/03-aging/color.majorclass.rds")

## all ##
pdf(paste0(outdir,"/integrated_with_class.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "class_label",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf(paste0(outdir,"/integrated_with_class.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "class_label",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color Joint-Cabernet samples only ##
merged.seuratobj.sct$class_label_color_our_only = merged.seuratobj.sct$class_label
merged.seuratobj.sct$class_label_color_our_only[which(merged.seuratobj.sct$source == "zeng")] = "zeng"
Idents(merged.seuratobj.sct) = merged.seuratobj.sct$class_label_color_our_only
pdf(paste0(outdir,"/integrated_with_class_color_our_only.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "class_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()
pdf(paste0(outdir,"/integrated_with_class_color_our_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "class_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

## color Zeng 10X samples only ##
merged.seuratobj.sct$class_label_color_zeng_only = merged.seuratobj.sct$class_label
merged.seuratobj.sct$class_label_color_zeng_only[which(merged.seuratobj.sct$source == "our")] = "our"
Idents(merged.seuratobj.sct) = merged.seuratobj.sct$class_label_color_zeng_only
pdf(paste0(outdir,"/integrated_with_class_color_zeng_only.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "class_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf(paste0(outdir,"/integrated_with_class_color_zeng_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "class_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## subclass ##
#inte.col = readRDS("/share/analysisdata/Methyl/public/TSO/yangfa/analysis/RNA_latest/01-Cluster_analysis/output/01-subclass_new.col_latest.rds")
#inte.col = c(inte.col,c("CA2-FC-IG Glut"="#89C75F","L6b CTX Glut"="#0C727C","Lamp5 Lhx6 Gaba"="#90D5E4","PAL-STR Gaba-Chol" ="#00ae9d","DG-PIR Ex IMN"="#1d953f","Vip Gaba"="#009ad6","L2/3 IT PIR-ENTl Glut"="#6E4B9E","L6 IT CTX Glut"="#AA0DFE","HPF CR Glut"='#e74c3c',"Pvalb chandelier Gaba"="#A6BDD7","STR D1 Sema5a Gaba" ="#B32851","OB-mi Frmd7 Gaba"='#5AC2F1FF',"OB Trdn Gaba"="#e4c6d0","OB Meis2 Thsd7b Gaba"="#f9906f","VLMC NN"="#ffc773","Sst Chodl Gaba"="#88c4e8","STR Prox1 Lhx6 Gaba"="#eb7f54","OT D3 Folh1 Gaba"="#815463","ABC NN" ="#253494","BAM NN" ="#FFFF00","Endo NN"="#d6ecf0","Peri NN"="#DEA0FD","Lymphoid NN"="#808080","SMC NN"="#bce672",'zeng'= 'lightgrey',"our"='lightgrey',"LA-BLA-BMA-PA Glut"="#666600","IT AON-TT-DP Glut"="#6A59EE"))
inte.col = readRDS("../../input/03-aging/color.subclass.rds")
## all ##
pdf(paste0(outdir,"/integrated_with_subclass.umap.pdf"),width = 20,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_label",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf(paste0(outdir,"/integrated_with_subclass.tsne.pdf"),width = 20,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass_label",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color Joint-Cabernet samples only ##
merged.seuratobj.sct$subclass_label_color_our_only = merged.seuratobj.sct$subclass_label
merged.seuratobj.sct$subclass_label_color_our_only[which(merged.seuratobj.sct$source == "zeng")] = "zeng"
Idents(merged.seuratobj.sct) = merged.seuratobj.sct$subclass_label_color_our_only
pdf(paste0(outdir,"/integrated_with_subclass_color_our_only.umap.pdf"),width = 20,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()
pdf(paste0(outdir,"/integrated_with_subclass_color_our_only.tsne.pdf"),width = 20,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

## color Zeng 10X samples only ##
merged.seuratobj.sct$subclass_label_color_zeng_only = merged.seuratobj.sct$subclass_label
merged.seuratobj.sct$subclass_label_color_zeng_only[which(merged.seuratobj.sct$source == "our")] = "our"
Idents(merged.seuratobj.sct) = merged.seuratobj.sct$subclass_label_color_zeng_only
pdf(paste0(outdir,"/integrated_with_subclass_color_zeng_only.umap.pdf"),width = 20,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf(paste0(outdir,"/integrated_with_subclass_color_zeng_only.tsne.pdf"),width = 20,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass_label_color_our_only",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## seurat cluster ##
pdf(paste0(outdir,"/seurat_cluster.integrated.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "seurat_clusters",label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()
pdf(paste0(outdir,"/seurat_cluster.integrated.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "seurat_clusters",label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

## source ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$source
levels(merged.seuratobj.sct)

inte.col = c("zeng"='#3498db',"our"='#e74c3c')
pdf(paste0(outdir,"/source.integrated.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "source",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()
pdf(paste0(outdir,"/source.integrated.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "source",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

## age ##
merged.seuratobj.sct$age[which(is.na(merged.seuratobj.sct$age))]='zeng'
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$age
levels(merged.seuratobj.sct)

inte.col = c("zeng"='#3498db',"old"='#e74c3c','young'='#FEE500')
pdf(paste0(outdir,"/age.integrated.umap.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "age",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()
pdf(paste0(outdir,"/age.integrated.tsne.pdf"),width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "age",cols = inte.col,label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

result_df <- data.frame(seurat_cluster = character(),major_class = character(), subclass = character(), source = character(), cell_number = numeric(), stringsAsFactors = FALSE)
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct@meta.data$seurat_clusters
for (ident_value in levels(merged.seuratobj.sct)) {
  cell_indices <- WhichCells(merged.seuratobj.sct, idents = ident_value)
  major_class <- unique(merged.seuratobj.sct@meta.data[cell_indices,"class_label"])
  subclass <- unique(merged.seuratobj.sct@meta.data[cell_indices,"subclass_label"])
  for(source in c("zeng","young","old")){
     cell_number = length(intersect(which(merged.seuratobj.sct$seurat_clusters == ident_value),which(merged.seuratobj.sct$age==source)))
     result_df <- rbind(result_df, data.frame(seurat_cluster = ident_value,major_class = major_class, subclass = subclass, source = source, cell_number = cell_number))
  }
  
}
print(result_df)

## save result ##
saveRDS(merged.seuratobj.sct,file=paste0(outdir,"/merged.seuratobj.sct.rds"))
saveRDS(our.seuratobj,file=paste0(outdir,"/our.seuratobj.min_dist_0.3.rds"))

##### 09. plot result of Joint-Cabernet samples #####
our.seuratobj <- RunUMAP(our.seuratobj, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=80)
## seurat cluster ##
pdf(paste0(outdir,"/seurat_cluster.integrated_our_cells_only.umap.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "seurat_clusters",label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()
pdf(paste0(outdir,"/seurat_cluster.integrated_our_cells_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "seurat_clusters",label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

## major class ##
inte.col = readRDS("../../input/03-aging/color.majorclass.rds")
## all ##
Idents(our.seuratobj) = our.seuratobj$lt_twice_class
pdf(paste0(outdir,"/major_class.integrated_our_cells_only.umap.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "lt_twice_class",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/major_class.integrated_our_cells_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "lt_twice_class",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

## color Joint-Cabernet young samples only ##
our.seuratobj$class_label_color_young_only = our.seuratobj$lt_twice_class
our.seuratobj$class_label_color_young_only[match(rownames(merged.seuratobj.sct@meta.data)[which(merged.seuratobj.sct$age == "old")],rownames(our.seuratobj@meta.data))] = "old"
Idents(our.seuratobj) = our.seuratobj$class_label_color_young_only
pdf(paste0(outdir,"/major_class.color_young_only.integrated_our_cells_only.umap.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "class_label_color_young_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(our.seuratobj))
dev.off()
pdf(paste0(outdir,"/major_class.color_young_only.integrated_our_cells_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "class_label_color_young_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(our.seuratobj))
dev.off()

## color Joint-Cabernet old samples only ##
our.seuratobj$class_label_color_old_only = our.seuratobj$lt_twice_class
our.seuratobj$class_label_color_old_only[match(rownames(merged.seuratobj.sct@meta.data)[which(merged.seuratobj.sct$age == "young")],rownames(our.seuratobj@meta.data))] = "young"
Idents(our.seuratobj) = our.seuratobj$class_label_color_old_only
#setdiff(unique(merged.seuratobj.sct$class_label),names(inte.col))
inte.col = c(inte.col,c("CNU-MGE GABA"="#FEE500","Vascular"='#00a6ac','old'= 'lightgrey',"young"='lightgrey'))
pdf(paste0(outdir,"/major_class.color_old_only.integrated_our_cells_only.umap.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "class_label_color_old_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/major_class.color_old_only.integrated_our_cells_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "class_label_color_old_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

## subclass ##
inte.col = readRDS("../../input/03-aging/color.subclass.rds")
Idents(our.seuratobj) = our.seuratobj$lt_twice_subclass
pdf(paste0(outdir,"/subclass.integrated_our_cells_only.umap.pdf"),width = 18,height = 13)
DimPlot(our.seuratobj, group.by = "lt_twice_subclass",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/subclass.integrated_our_cells_only.tsne.pdf"),width = 18,height = 13)
TSNEPlot(our.seuratobj, group.by = "lt_twice_subclass",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

## color Joint-Cabernet young samples only ##
our.seuratobj$subclass_label_color_young_only = our.seuratobj$lt_twice_subclass
our.seuratobj$subclass_label_color_young_only[match(rownames(merged.seuratobj.sct@meta.data)[which(merged.seuratobj.sct$age == "old")],rownames(our.seuratobj@meta.data))] = "old"
Idents(our.seuratobj) = our.seuratobj$subclass_label_color_young_only
pdf(paste0(outdir,"/subclass.color_young_only.integrated_our_cells_only.umap.pdf"),width = 18,height = 13)
DimPlot(our.seuratobj, group.by = "subclass_label_color_young_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(our.seuratobj))
dev.off()
pdf(paste0(outdir,"/subclass.color_young_only.integrated_our_cells_only.tsne.pdf"),width = 18,height = 13)
TSNEPlot(our.seuratobj, group.by = "subclass_label_color_young_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(our.seuratobj))
dev.off()

## color Joint-Cabernet old samples only ##
our.seuratobj$subclass_label_color_old_only = our.seuratobj$lt_twice_subclass
our.seuratobj$subclass_label_color_old_only[match(rownames(merged.seuratobj.sct@meta.data)[which(merged.seuratobj.sct$age == "young")],rownames(our.seuratobj@meta.data))] = "young"
Idents(our.seuratobj) = our.seuratobj$subclass_label_color_old_only
pdf(paste0(outdir,"/subclass.color_old_only.integrated_our_cells_only.umap.pdf"),width = 18,height = 13)
DimPlot(our.seuratobj, group.by = "subclass_label_color_old_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/subclass.color_old_only.integrated_our_cells_only.tsne.pdf"),width = 18,height = 13)
TSNEPlot(our.seuratobj, group.by = "subclass_label_color_old_only",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

## age ##
our.seuratobj$age=merged.seuratobj.sct@meta.data[rownames(our.seuratobj@meta.data),"age"]
Idents(our.seuratobj) <- our.seuratobj$age
levels(our.seuratobj)

## all ##
inte.col = c('young'='#3498db',"old"='#e74c3c')
pdf(paste0(outdir,"/age.integrated_our_cells_only.umap.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "age",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()
pdf(paste0(outdir,"/age.integrated_our_cells_only.tsne.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "age",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

Idents(our.seuratobj) = our.seuratobj$age
inte.col=c('young'='#3498db',"old"='#9e9e9e')
pdf(paste0(outdir,"/our_age_color_young_only_distribution.min_dist_0.5.alpha_0.5.pt.size_1.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "age",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,alpha = 0.5,order=c("young","old"))
dev.off()

inte.col=c('young'='#9e9e9e',"old"='#e74c3c')
pdf(paste0(outdir,"/our_age_color_old_only_distribution.min_dist_0.5.alpha_0.5.pt.size_1.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "age",cols = inte.col,label=TRUE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,alpha = 0.5)
dev.off()

## without label ##
inte.col=c('young'='#3498db',"old"='#9e9e9e')
pdf(paste0(outdir,"/our_age_color_young_only_distribution.min_dist_0.5.alpha_0.5.pt.size_1.without_label.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "age",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,alpha = 0.5,order=c("young","old"))
dev.off()

inte.col=c('young'='#9e9e9e',"old"='#e74c3c')
pdf(paste0(outdir,"/our_age_color_old_only_distribution.min_dist_0.5.alpha_0.5.pt.size_1.without_label.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "age",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,alpha = 0.5)
dev.off()

## subclass ##
inte.col = readRDS("../../input/02-aging/color.subclass.rds")
Idents(our.seuratobj) = our.seuratobj$lt_twice_subclass
pdf(paste0(outdir,"/subclass.integrated_our_cells_only.umap.without_label.pdf"),width = 18,height = 13)
DimPlot(our.seuratobj, group.by = "lt_twice_subclass",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/subclass.integrated_our_cells_only.tsne.without_label.pdf"),width = 18,height = 13)
TSNEPlot(our.seuratobj, group.by = "lt_twice_subclass",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

## major class ##
inte.col = readRDS("../../input/03-aging/color.majorclass.rds")
Idents(our.seuratobj) = our.seuratobj$lt_twice_class
pdf(paste0(outdir,"/major_class.integrated_our_cells_only.umap.without_label.pdf"),width = 15,height = 13)
DimPlot(our.seuratobj, group.by = "lt_twice_class",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()
pdf(paste0(outdir,"/major_class.integrated_our_cells_only.tsne.without_label.pdf"),width = 15,height = 13)
TSNEPlot(our.seuratobj, group.by = "lt_twice_class",cols = inte.col,label=FALSE, pt.size =1,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(our.seuratobj)))
dev.off()

##### subclass sample file list generate 
# paired_sample = read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/old_mouse_sample_info/01_Sample_info/RNA_DNA_match_name_QC.v20241009.csv")
paired_sample = read.csv("../../input/03-aging/RNA_DNA_match_name_QC_older_than_P70.csv")
metainfo = our.seuratobj@meta.data
rownames(metainfo)[which(!is.na(str_match(rownames(metainfo),"@@")))] = unlist(lapply(rownames(metainfo)[which(!is.na(str_match(rownames(metainfo),"@@")))], function(x) strsplit(x,"@@_")[[1]][2]))
#young_match.df = read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/sample_info/01_Sampleinfo/RNA_DNA_match_name_QC_44608.csv",header=T)
young_match.df = read.csv("../../input/01-youth/RNA_DNA_match_name_QC_class_label.csv",header=T)

## generate allc path for all subclassed all ages ##
# setwd("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/01.merge_allc/subclass_sample_file_list")
for(cl in unique(metainfo$lt_twice_subclass)){
    for(year in c("old","young")){
        if(year == "old"){
           path = ourdir
           for(type in c("5hmC","5mC")){
            select.df = subset(metainfo,lt_twice_subclass== cl & age == year)
            if(type == "5hmC"){
                sample_list = paired_sample[intersect(match(rownames(select.df),paired_sample$RNA),which(paired_sample$total_QC == 1)),"hmC"]
            }else{
                sample_list = paired_sample[intersect(match(rownames(select.df),paired_sample$RNA),which(paired_sample$total_QC == 1)),"mC"]
            }
            #write.table(paste0(path,type,"/allc_",sample_list,".mm10.dna.tsv.gz"),file=paste0("subclass.",str_replace_all(str_replace_all(cl, " ", "."), "/", "."),"_",year,"_",type,".txt"),quote=F,row.names=F,col.names=F)
            write.table(paste0(path,type,"/allc_",sample_list,".mm10.dna.tsv.gz"),file=paste0("../../output/02-aging/01-merge_allc/subclass.",str_replace_all(str_replace_all(cl, " ", "."), "/", "."),"_",year,"_",type,".txt"),quote=F,row.names=F,col.names=F)
           }
        }else{
          path = ourdir
          for(type in c("5hmC","5mC")){
            select.df = subset(metainfo,lt_twice_subclass== cl & age == year)
            if(type == "5hmC"){
                sample_list = young_match.df[intersect(match(rownames(select.df),young_match.df$RNA),which(young_match.df$total_QC == 1)),"hmC"]
            }else{
                sample_list = young_match.df[intersect(match(rownames(select.df),young_match.df$RNA),which(young_match.df$total_QC == 1)),"mC"]
            }
            #write.table(paste0(path,type,"/allc_",sample_list,".mm10.dna.tsv.gz"),file=paste0("subclass.",str_replace_all(str_replace_all(cl, " ", "."), "/", "."),"_",year,"_",type,".txt"),quote=F,row.names=F,col.names=F)
            write.table(paste0(path,type,"/allc_",sample_list,".mm10.dna.tsv.gz"),file=paste0("../../output/02-aging/01-merge_allc/subclass.",str_replace_all(str_replace_all(cl, " ", "."), "/", "."),"_",year,"_",type,".txt"),quote=F,row.names=F,col.names=F)
           }
        }
        
    }
}

## check ##
rownames(paired_sample) = paired_sample$RNA
paired_sample = paired_sample[rownames(metainfo),]
metainfo$RNA_QC = paired_sample$RNA_QC
metainfo$total_QC = paired_sample$total_QC
metainfo$mC_QC = paired_sample$mC_QC
dim(metainfo[metainfo$age == "old" & metainfo$lt_twice_subclass == "Astro-TE NN" & metainfo$total_QC == 1,])
# 1718

## save results ##
saveRDS(our.seuratobj,file=paste0(outdir,"/our.seuratobj.min_dist_0.5.rds"))
save(our.seuratobj,merged.seuratobj.sct,file=paste0(outdir,"/seurat_obj.our_integrated.RData"))

## generate new metainfo ##
info=read.csv("../../input/03-aging/total_info.csv",header=T)
rownames(info) = info$SampleID
info = info[rownames(metainfo),]

metainfo$Batch = info$Batch
metainfo$mC = info$mC
metainfo$hmC = info$hmC
metainfo$RNA_QC = info$RNA_QC
metainfo$mC_QC = info$mC_QC
metainfo$hmC_QC = info$hmC_QC
metainfo$total_QC = info$total_QC

metainfo$three_class = NA
metainfo$three_class[which(!is.na(str_match(metainfo$lt_twice_subclass,"Glut")))] = "Glut"
metainfo$three_class[which(!is.na(str_match(metainfo$lt_twice_subclass,"Gaba")))] = "Gaba"
metainfo$three_class[which(!is.na(str_match(metainfo$lt_twice_subclass,"NN")))] = "NN"

saveRDS(metainfo,file="../../output/03-aging/metainfo.250115.rds")
# write.csv(metainfo,file="mouse_young_and_old.metainfo.csv",quote=F,row.names=T,col.names=T)
write.csv(metainfo,file="../../output/03-aging/mouse_young_and_old.metainfo.csv",quote=F,row.names=T,col.names=T)