##### 01. import packages #####
library(getopt)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

##### 02. parameter passing #####
arguments = matrix(c(
  'help', 'h', 0, "logical", "",
  'neuron_type', 'n', 1, "character", "",
  'label', 'l', 1, "character", ""
), byrow=TRUE, ncol=5)
args = getopt(arguments)

if (!is.null(args$help) || is.null(args$neuron_type) || is.null(args$label)) {
  cat(paste(getopt(arguments, usage = T), "\n"))
  q()
}

##### 03. set working directory #####
setwd("./")

##### 04. mapping #####
#merged.seuratobj.sct <- readRDS(paste0(args$neuron_type,"_obj.label_transfer_",args$label,".rds"))
merged.seuratobj.sct <- readRDS(paste0("../output/02.Young_Mouse.Brain_slice/integration/",args$neuron_type,"_obj.label_transfer_",args$label,".rds"))

## extract UMAP loci ##
TSNE_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  # extract final UMAP loci

## calculate the euclidean distance ##
dist_matrix <- as.dist(dist(TSNE_RNA))
dist_matrix_mat <- as.matrix(dist_matrix)  

dist_matrix_mat <- dist_matrix_mat[which(!is.na(str_match(rownames(dist_matrix_mat),"Mouses"))),]   # reserve Joint-Cabernet RNA sample as row
dist_matrix_mat <- dist_matrix_mat[,which(is.na(str_match(colnames(dist_matrix_mat),"Mouses")))]   # reserve MERFISH sample as col

metainfo <- merged.seuratobj.sct@meta.data
metainfo_subset <- metainfo[which(is.na(str_match(rownames(metainfo),"Mouses"))),]

## calculate the nearst Zhuang MERFISH cell for each Joint-Cabernet slice cell ##
result_df <- data.frame(idents = character(), nearst = character(), stringsAsFactors = FALSE)
for(i in 1:dim(dist_matrix_mat)[1]){
    cl <- metainfo[rownames(dist_matrix_mat)[i],"subclass"] # extract subclass info of Joint-Cabernet slice cell
    region <- metainfo[rownames(dist_matrix_mat)[i],"major_brain_region_v2"] # extract brain_region info of Joint-Cabernet slice cell
    cells <- rownames(metainfo[intersect(which(metainfo$subclass == cl),which(metainfo$major_brain_region_v2 == region)),]) # limite the candidate cells to the same cell type in the same brain region.
    temp <- dist_matrix_mat[i,match(cells,colnames(dist_matrix_mat))] # extract distance
    indices <- order(temp,decreasing = FALSE)  # extract the nearst sample index
    result_df <- rbind(result_df, data.frame(idents = rownames(dist_matrix_mat)[i], nearst = paste(names(temp)[indices[1:1]],collapse=";"))) # record
}

## extract spatial loci of Zhuang MERFISH cells nearst to each Joint-Cabernet slice cell ##
## Zhuang MERFISH metainfo
loc.df <- read.csv("../04.data/02.metainfo/02.Young_Mouse.Brain_slice/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
rownames(loc.df) <- loc.df$cell_label

x.v <- c()
y.v <- c()
for(i in 1:dim(result_df)[1]){
    cells <- unlist(lapply(result_df[i,"nearst"], function(x) strsplit(x,";")[[1]]))
    x.v <- c(x.v,mean(loc.df[cells,"x"]))
    y.v <- c(y.v,mean(loc.df[cells,"y"]))
}

result_df$x <- x.v
result_df$y <- y.v 


order.v <- rownames(merged.seuratobj.sct@meta.data)
colnames(TSNE_RNA) <- c("x","y")   # change to "umap_1","umap_2" later
TSNE_RNA$x <- NA
TSNE_RNA$y <- NA

## Zhuang MERFISH cell spatial loci substitution ##
TSNE_RNA$x[which(is.na(str_match(colnames(dist_matrix_mat),"Mouses")))] <- loc.df[rownames(TSNE_RNA)[which(is.na(str_match(colnames(dist_matrix_mat),"Mouses")))],'x']
TSNE_RNA$y[which(is.na(str_match(colnames(dist_matrix_mat),"Mouses")))] <- loc.df[rownames(TSNE_RNA)[which(is.na(str_match(colnames(dist_matrix_mat),"Mouses")))],'y']

## Joint-Cabernet slice cell spatial loci substitution ##
rownames(result_df) <- result_df$ide
TSNE_RNA$x[which(!is.na(str_match(rownames(TSNE_RNA),"Mouses")))] <- result_df[rownames(TSNE_RNA)[which(!is.na(str_match(rownames(TSNE_RNA),"Mouses")))],'x']
TSNE_RNA$y[which(!is.na(str_match(rownames(TSNE_RNA),"Mouses")))] <- result_df[rownames(TSNE_RNA)[which(!is.na(str_match(rownames(TSNE_RNA),"Mouses")))],'y']

TSNE_RNA <- TSNE_RNA[order.v,] # order

colnames(TSNE_RNA) <- c("umap_1","umap_2")
TSNE_RNA <- as.matrix(TSNE_RNA)
## replace to umap dimension ##
merged.seuratobj.sct@reductions$umap@cell.embeddings <- TSNE_RNA # assign to UMAP for plot

## save result ##
saveRDS(merged.seuratobj.sct,file=paste0("../output/02.Young_Mouse.Brain_slice/map/",args$neuron_type,"_obj.label_transfer_",args$label,".loci_transfer.the_nearst_1_cell.rds"))

