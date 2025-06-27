##### 00. enviroment #####
#conda activate /share/analysisdata/Methyl/public/rna

##### 01.import packages #####
library(Seurat)
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library"))
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)

##### 02.read data and subset data whose Z coordinates and brain region consistent with Joint-Cabernet slice data #####
## 02.01 subset MERFISH data of Zhuang by z with threshold of 7.33 ##

## rds ##
zhuang.obj <- readRDS("../04.data/03.download_data/Zhuang_raw_counts.rds")
# downloaded from https://datasets.cellxgene.cziscience.com/f9f6cc29-7d06-47be-a798-e4e3b36a86b2.rds

## annotation ##
anno.df <- read.csv("../04.data/02.metainfo/02.Young_Mouse.Brain_slice/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
idx.v <- intersect(which(anno.df$z > 7.33),which(anno.df$z < 7.34))
select.v <-  loc[idx.v,'cell_label']

zhuang.obj$select_v <- "N"
zhuang.obj$select_v[which(rownames(zhuang.obj@meta.data) %in% select.v)] <- "Y"
seurat_subset <- subset(zhuang.obj, select_v == "Y")

## 02.02 subset MERFISH data of Zhuang by major_brain_region ###
seurat_subset@meta.data$select <- "N"
seurat_subset@meta.data$select[c(which(seurat_subset@meta.data$major_brain_region == "Isocortex"),which(seurat_subset@meta.data$major_brain_region == "Hippocampus"))] <- "Y"
seurat_subset_region <- subset(seurat_subset, select == "Y")

##### 03.save rds for further analysis #####
saveRDS(seurat_subset_region,file="../output/02.Young_Mouse.Brain_slice/subset.z_axis_located_on_7.33.cortex_and_hippo.rds")
