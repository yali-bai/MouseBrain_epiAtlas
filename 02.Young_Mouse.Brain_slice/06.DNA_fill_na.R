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


merged.seuratobj.sct <- readRDS("merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
gene.v <- rownames(merged.seuratobj.sct@assays$RNA$counts)
metainfo = merged.seuratobj.sct@meta.data

paired_sampleinfo = read.csv("../03.data/02.metainfo/02.Young_Mouse.Brain_slice/RNA_DNA_match_name_QC_class_label_young.brain_slice.add_celltype.csv",header=T)

meta_info.df <- merged.seuratobj.sct@meta.data
meta_info.df <- meta_info.df[which(!is.na(str_match(rownames(meta_info.df),"TSO"))),]
fill_na <- function(matrix,type){
    #matrix <- promoter_hmCG_mCG
    matrix_fill <- matrix
    #type <- "5mC"
    if(type=="5mC"){
        coln = "mC_SampleID"
    }else{
        coln = "hmC_SampleID"
    }
    rownames(paired_sampleinfo) = paired_sampleinfo$RNA_SampleID
    for(i in 1:dim(matrix)[1]){
        na.idx <- which(is.na(matrix[i,]))
        for(cl in unique(meta_info.df$subclass)){
            methy_v <- mean(as.numeric(matrix[i,which(colnames(matrix) %in% paired_sampleinfo[rownames(meta_info.df[which(meta_info.df$subclass == cl),]),coln])]),na.rm=TRUE)
            matrix_fill[i,intersect(na.idx,match(paired_sampleinfo[rownames(meta_info.df[which(meta_info.df$subclass == cl),]),coln],colnames(matrix)))] <- methy_v
        }
    }
    return(matrix_fill)
}

tmp_DNA <- MuDataSeurat::ReadH5AD("5hmC_genebody.CG.pass_total_QC.h5ad")
gene_hmCG <- as.data.frame(tmp_DNA@assays$RNA@counts) 
rownames(gene_hmCG) <- unlist(lapply(rownames(gene_hmCG), function(x) strsplit(x,"\\.")[[1]][1]))
gene_hmCG <- gene_hmCG[match(intersect(gene.v,rownames(gene_hmCG)),rownames(gene_hmCG)),]

tmp_DNA <- MuDataSeurat::ReadH5AD("5mC_genebody.CG.pass_total_QC.h5ad")
gene_hmCG_mCG <- as.data.frame(tmp_DNA@assays$RNA@counts) 
rownames(gene_hmCG_mCG) <- unlist(lapply(rownames(gene_hmCG_mCG), function(x) strsplit(x,"\\.")[[1]][1]))
gene_hmCG_mCG <- gene_hmCG_mCG[match(gene.v,rownames(gene_hmCG_mCG)),]

gene_hmCG_mCG = gene_hmCG_mCG[rownames(gene_hmCG),]
rownames(paired_sampleinfo) = paired_sampleinfo$hmC_SampleID
colnames(gene_hmCG_mCG) = str_replace_all(colnames(gene_hmCG_mCG),"allc_","")
colnames(gene_hmCG) = str_replace_all(colnames(gene_hmCG),"allc_","")
gene_hmCG_mCG = gene_hmCG_mCG[,paired_sampleinfo[colnames(gene_hmCG),"mC_SampleID"]]
gene_mCG = gene_hmCG_mCG - gene_hmCG

retain_mC = intersect(colnames(gene_mCG),paired_sampleinfo[paired_sampleinfo$total_QC == 1,"mC_SampleID"])
retain_hmC = intersect(colnames(gene_hmCG),paired_sampleinfo[paired_sampleinfo$total_QC == 1,"hmC_SampleID"])

gene_mCG = gene_mCG[,retain_mC]
gene_hmCG_mCG = gene_hmCG_mCG[,retain_mC]
gene_hmCG = gene_hmCG[,retain_hmC]

gene_hmCG_fill <- fill_na(gene_hmCG,"5hmC")
gene_hmCG_mCG_fill <- fill_na(gene_hmCG_mCG,"5mC")
gene_mCG_fill <- fill_na(gene_mCG,"5mC")

save(gene_hmCG_fill,gene_hmCG_mCG_fill,gene_mCG_fill,file="DNA_fill_na.filter_by_total_QC.RData")









