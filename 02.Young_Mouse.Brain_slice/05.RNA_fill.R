## packages
library(Seurat)
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(MuDataSeurat)
library(grid)
library(cowplot)
library(data.table) 

## our seurat object
merged.seuratobj.sct <- readRDS("merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")

fill_na_RNA <- function(matrix,type){
    matrix_fill <- matrix
    for(i in 1:dim(matrix)[1]){
        na.idx <- which(matrix[i,] == 0)
        for(cl in unique(merged.seuratobj.sct@meta.data$subclass)){
            methy_v <- mean(as.numeric(matrix[i,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == cl),])]),na.rm=TRUE)
            matrix_fill[i,intersect(na.idx,match(rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == cl),]),colnames(matrix)))] <- methy_v
        }
    }
    return(matrix_fill)
}
RNA_fill <- fill_na_RNA(as.data.frame(merged.seuratobj.sct@assays$RNA$data))

saveRDS(RNA_fill,file = "RNA_fill_na.rds")