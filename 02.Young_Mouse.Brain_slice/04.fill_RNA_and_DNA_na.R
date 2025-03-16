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

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02. change working path #####
# setwd(paste0(indir,"")

## our seurat object ##
merged.seuratobj.sct <- readRDS("../output/02-slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
gene.v <- rownames(merged.seuratobj.sct@assays$RNA$counts)
metainfo = merged.seuratobj.sct@meta.data

## sample paired info ##
#paired_sampleinfo = read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/sample_info/01_Sampleinfo/RNA_DNA_match_name_QC_MERFISH.csv",header=T)
paired_sampleinfo = read.csv("../input/02-slice/RNA_DNA_match_name_QC_MERFISH.csv",header=T)

## load methyl matrix ##
## genebody ##
## hmC_CGN ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5hmC_genebody.CG.pass_hmC_QC.h5ad"))
geneslop2k_hmC_CGN <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5hmC_CGN.unnormalized.geneslop2k.csv.gz",row.names=1,header=T)
rownames(geneslop2k_hmC_CGN) <- unlist(lapply(rownames(geneslop2k_hmC_CGN), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(geneslop2k_hmC_CGN) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
geneslop2k_hmC_CGN <- geneslop2k_hmC_CGN[match(intersect(gene.v,rownames(geneslop2k_hmC_CGN)),rownames(geneslop2k_hmC_CGN)),]

## hmC_CHN ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5hmC_genebody.CH.pass_hmC_QC.h5ad"))
geneslop2k_hmC_CHN <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5hmC_CHN.unnormalized.geneslop2k.csv.gz",row.names=1,header=T)
rownames(geneslop2k_hmC_CHN) <- unlist(lapply(rownames(geneslop2k_hmC_CHN), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(geneslop2k_hmC_CHN) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
geneslop2k_hmC_CHN <- geneslop2k_hmC_CHN[match(gene.v,rownames(geneslop2k_hmC_CHN)),]

## hmCG_mCG ##
geneslop2k_hmCG_mCG <- fread(paste0(indir,"/5mC_CGN.genebody.csv.gz"),header=TRUE,sep=",",data.table = FALSE)
rownames(geneslop2k_hmCG_mCG) = geneslop2k_hmCG_mCG$V1
geneslop2k_hmCG_mCG = geneslop2k_hmCG_mCG[,-which(colnames(geneslop2k_hmCG_mCG) == "V1")]
rownames(geneslop2k_hmCG_mCG) <- unlist(lapply(rownames(geneslop2k_hmCG_mCG), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(geneslop2k_hmCG_mCG) <- unlist(lapply(colnames(geneslop2k_hmCG_mCG), function(x) strsplit(x,"allc_")[[1]][2]))
geneslop2k_hmCG_mCG <- geneslop2k_hmCG_mCG[match(gene.v,rownames(geneslop2k_hmCG_mCG)),]

## hmCH_mCH ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5mC_genebody.CH.pass_mC_QC.h5ad"))
geneslop2k_hmCH_mCH <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5mC_CHN.unnormalized.geneslop2k.csv.gz",row.names=1,header=T)
rownames(geneslop2k_hmCH_mCH) <- unlist(lapply(rownames(geneslop2k_hmCH_mCH), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(geneslop2k_hmCH_mCH) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
geneslop2k_hmCH_mCH <- geneslop2k_hmCH_mCH[match(gene.v,rownames(geneslop2k_hmCH_mCH)),]

## mC_CGN (mC = hmC_mC - hmC) ##
hmC.idx = match(colnames(geneslop2k_hmC_CGN),paired_sampleinfo$hmC)
mC.idx = match(colnames(geneslop2k_hmCG_mCG),paired_sampleinfo$mC)
intersect.idx = intersect(hmC.idx,mC.idx)
intersect.idx = intersect(intersect.idx,which(paired_sampleinfo$total_QC == 1))
geneslop2k_hmC_CGN <- geneslop2k_hmC_CGN[rownames(geneslop2k_hmCG_mCG),]
geneslop2k_hmC_CGN <- geneslop2k_hmC_CGN[,paired_sampleinfo$hmC[intersect.idx]]
geneslop2k_hmCG_mCG <- geneslop2k_hmCG_mCG[,paired_sampleinfo$mC[intersect.idx]]
geneslop2k_mC_CGN  <- geneslop2k_hmCG_mCG - geneslop2k_hmC_CGN
geneslop2k_mC_CGN <- geneslop2k_mC_CGN[match(gene.v,rownames(geneslop2k_mC_CGN)),]

## mC_CHN ##
hmC.idx = match(colnames(geneslop2k_hmC_CHN),paired_sampleinfo$hmC)
mC.idx = match(colnames(geneslop2k_hmCH_mCH),paired_sampleinfo$mC)
intersect.idx = intersect(hmC.idx,mC.idx)
intersect.idx = intersect(intersect.idx,which(paired_sampleinfo$total_QC == 1))
geneslop2k_hmC_CHN <- geneslop2k_hmC_CHN[rownames(geneslop2k_hmCH_mCH),]
geneslop2k_hmC_CHN <- geneslop2k_hmC_CHN[,paired_sampleinfo$hmC[intersect.idx]]
geneslop2k_hmCH_mCH <- geneslop2k_hmCH_mCH[,paired_sampleinfo$mC[intersect.idx]]
geneslop2k_mC_CHN <- geneslop2k_hmCH_mCH - geneslop2k_hmC_CHN
geneslop2k_mC_CHN <- geneslop2k_mC_CHN[match(gene.v,rownames(geneslop2k_mC_CHN)),]

## promoter ##
## hmC_CGN ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5hmC_promoter.CG.pass_hmC_QC.h5ad"))
promoter_hmC_CGN <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5hmC_CGN.unnormalized.promoter.csv.gz",row.names=1,header=T)
rownames(promoter_hmC_CGN) <- unlist(lapply(rownames(promoter_hmC_CGN), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(promoter_hmC_CGN) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
promoter_hmC_CGN <- promoter_hmC_CGN[match(gene.v,rownames(promoter_hmC_CGN)),]

## hmC_CHN ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5hmC_promoter.CH.pass_hmC_QC.h5ad"))
promoter_hmC_CHN <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5hmC_CHN.unnormalized.promoter.csv.gz",row.names=1,header=T)
rownames(promoter_hmC_CHN) <- unlist(lapply(rownames(promoter_hmC_CHN), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(promoter_hmC_CHN) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
promoter_hmC_CHN <- promoter_hmC_CHN[match(gene.v,rownames(promoter_hmC_CHN)),]

## hmCG_mCG ##
promoter_hmCG_mCG <- fread(paste0(indir,"/5mC_CGN.promoter.csv.gz"),header=TRUE,sep=",",data.table = FALSE)
rownames(promoter_hmCG_mCG) = promoter_hmCG_mCG$V1
promoter_hmCG_mCG = promoter_hmCG_mCG[,-which(colnames(promoter_hmCG_mCG) == "V1")]
rownames(promoter_hmCG_mCG) <- unlist(lapply(rownames(promoter_hmCG_mCG), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(promoter_hmCG_mCG) <- unlist(lapply(colnames(promoter_hmCG_mCG), function(x) strsplit(x,"allc_")[[1]][2]))
promoter_hmCG_mCG <- promoter_hmCG_mCG[match(gene.v,rownames(promoter_hmCG_mCG)),]

## hmCH_mCH ##
tmp_DNA <- MuDataSeurat::ReadH5AD(paste0(indir,"/5mC_promoter.CH.pass_mC_QC.h5ad"))
promoter_hmCH_mCH <- as.data.frame(tmp_DNA@assays$RNA@counts) #read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/DNA/MERFISH/raw_fraction/5mC_CHN.unnormalized.promoter.csv.gz",row.names=1,header=T)
rownames(promoter_hmCH_mCH) <- unlist(lapply(rownames(promoter_hmCH_mCH), function(x) strsplit(x,"\\.")[[1]][1]))
colnames(promoter_hmCH_mCH) <- unlist(lapply(rownames(tmp_DNA@meta.data), function(x) strsplit(x,"allc_")[[1]][2]))
promoter_hmCH_mCH <- promoter_hmCH_mCH[match(gene.v,rownames(promoter_hmCH_mCH)),]

## mC_CGN ##
hmC.idx = match(colnames(promoter_hmC_CGN),paired_sampleinfo$hmC)
mC.idx = match(colnames(promoter_hmCG_mCG),paired_sampleinfo$mC)
intersect.idx = intersect(hmC.idx,mC.idx)
intersect.idx = intersect(intersect.idx,which(paired_sampleinfo$total_QC == 1))
promoter_hmC_CGN <- promoter_hmC_CGN[rownames(promoter_hmCG_mCG),]
promoter_hmC_CGN <- promoter_hmC_CGN[,paired_sampleinfo$hmC[intersect.idx]]
promoter_hmCG_mCG <- promoter_hmCG_mCG[,paired_sampleinfo$mC[intersect.idx]]
promoter_mC_CGN  <- promoter_hmCG_mCG - promoter_hmC_CGN
promoter_mC_CGN <- promoter_mC_CGN[match(gene.v,rownames(promoter_mC_CGN)),]

## mC_CHN ##
hmC.idx = match(colnames(promoter_hmC_CHN),paired_sampleinfo$hmC)
mC.idx = match(colnames(promoter_hmCH_mCH),paired_sampleinfo$mC)
intersect.idx = intersect(hmC.idx,mC.idx)
intersect.idx = intersect(intersect.idx,which(paired_sampleinfo$total_QC == 1))
promoter_hmC_CHN <- promoter_hmC_CHN[rownames(promoter_hmCH_mCH),]
promoter_hmC_CHN <- promoter_hmC_CHN[,paired_sampleinfo$hmC[intersect.idx]]
promoter_hmCH_mCH <- promoter_hmCH_mCH[,paired_sampleinfo$mC[intersect.idx]]
promoter_mC_CHN <- promoter_hmCH_mCH - promoter_hmC_CHN
promoter_mC_CHN <- promoter_mC_CHN[match(gene.v,rownames(promoter_mC_CHN)),]

## fill DNA NA ##
meta_info.df <- merged.seuratobj.sct@meta.data
meta_info.df <- meta_info.df[which(!is.na(str_match(rownames(meta_info.df),"Mouse"))),]
## define function which can fill na values as subclass mean ##
fill_na <- function(matrix,type){
    #matrix <- promoter_hmCG_mCG
    #type <- "5mC"
    
    matrix_fill <- matrix
    if(type=="5mC"){
        coln = "mC"
    }else{
        coln = "hmC"
    }
    rownames(paired_sampleinfo) = paired_sampleinfo$RNA
    for(i in 1:dim(matrix)[1]){
        na.idx <- which(is.na(matrix[i,]))
        for(cl in unique(meta_info.df$subclass)){
            methy_v <- mean(as.numeric(matrix[i,which(colnames(matrix) %in% paired_sampleinfo[rownames(meta_info.df[which(meta_info.df$subclass == cl),]),coln])]),na.rm=TRUE)
            matrix_fill[i,intersect(na.idx,match(paired_sampleinfo[rownames(meta_info.df[which(meta_info.df$subclass == cl),]),coln],colnames(matrix)))] <- methy_v
        }
    }
    return(matrix_fill)
}

## run function ##
## genebody ## 
geneslop2k_hmC_CGN_fill <- fill_na(geneslop2k_hmC_CGN,"5hmC")
geneslop2k_hmC_CHN_fill <- fill_na(geneslop2k_hmC_CHN,"5hmC") 
geneslop2k_hmCG_mCG_fill <- fill_na(geneslop2k_hmCG_mCG,"5mC")
geneslop2k_hmCH_mCH_fill <- fill_na(geneslop2k_hmCH_mCH,"5mC")
geneslop2k_mC_CGN_fill <- fill_na(geneslop2k_mC_CGN,"5mC")
geneslop2k_mC_CHN_fill <- fill_na(geneslop2k_mC_CHN,"5mC")

## promoter ##
promoter_hmC_CGN_fill <- fill_na(promoter_hmC_CGN,"5hmC")
promoter_hmC_CHN_fill <- fill_na(promoter_hmC_CHN,"5hmC") 
promoter_hmCG_mCG_fill <- fill_na(promoter_hmCG_mCG,"5mC")
promoter_hmCH_mCH_fill <- fill_na(promoter_hmCH_mCH,"5mC")
promoter_mC_CGN_fill <- fill_na(promoter_mC_CGN,"5mC")
promoter_mC_CHN_fill <- fill_na(promoter_mC_CHN,"5mC")

## fill RNA na ##
## define function which can fill zero values as subclass mean ##
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

## run function ##
RNA_fill <- fill_na_RNA(as.data.frame(merged.seuratobj.sct@assays$RNA$data))

## save result ##
save(geneslop2k_hmC_CGN_fill,geneslop2k_hmC_CHN_fill,geneslop2k_hmCG_mCG_fill,geneslop2k_hmCH_mCH_fill,geneslop2k_mC_CGN_fill,geneslop2k_mC_CHN_fill,
    promoter_hmC_CGN_fill,promoter_hmC_CHN_fill,promoter_hmCG_mCG_fill,promoter_hmCH_mCH_fill,promoter_mC_CGN_fill,promoter_mC_CHN_fill,
    RNA_fill, file="../output/02-slice/RNA_DNA_fill_na.20240925.RData")


