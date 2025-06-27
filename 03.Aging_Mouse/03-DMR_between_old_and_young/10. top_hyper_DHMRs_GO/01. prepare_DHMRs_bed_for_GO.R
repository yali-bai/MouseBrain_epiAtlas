##### 01.import packages #####
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library","/share/analysisdata/Methyl/public/rna/lib/R/library"))
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(ggpointdensity) 
library(cowplot)
library(ggtext)
library(ggpubr)
library(ggunchained)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
# setwd(paste0(outdir,"")

##### 03.read subclass order file #####
subclass_order = readRDS("../../../04.data/04.config_files/order.subclass.rds")

##### 04.get top hyper DHMRs position #####
## read significant DHMRs file in ##
DHMR_result=readRDS("../../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")
head(DHMR_result)

## select hyper DHMRs ##
DHMR_result.sig = subset(DHMR_result, diff > 0.05 )

## sort by subclass and diff ##
DHMR_result.sig$cluster = factor(DHMR_result.sig$cluster,levels = intersect(subclass_order,unique(DHMR_result.sig$cluster)))
DHMR_result.sig.sorted <- DHMR_result.sig[order(DHMR_result.sig$cluster, -DHMR_result.sig$diff), ]
DHMR_result.sig.sorted[1:5,]

## extract chr info, the start and end info have been extracted when filtering length ##
DHMR_result.sig.sorted$chr = unlist(lapply(DHMR_result.sig.sorted$chrom, function(x) strsplit(x,'_')[[1]][1]))
unique(DHMR_result.sig.sorted$cluster)

## assign three class ##
DHMR_result.sig.sorted$three_class = "Non-neuron"
DHMR_result.sig.sorted$three_class[which(!is.na(str_match(DHMR_result.sig.sorted$cluster,"Glut")))] = "Exc"
DHMR_result.sig.sorted$three_class[which(!is.na(str_match(DHMR_result.sig.sorted$cluster,"Gaba")))] = "Inh"

## all hyper DHMRs position of all subclasses ##
if (!dir.exists(paste0(outdir,"/subclass/all"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/subclass/all"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$cluster)){
    df = data.frame(chr=DHMR_result.sig.sorted[DHMR_result.sig.sorted$cluster == cl,"chr"],start=DHMR_result.sig.sorted[DHMR_result.sig.sorted$cluster == cl,"start"],end=DHMR_result.sig.sorted[DHMR_result.sig.sorted$cluster == cl,"end"])
    write.table(df,file=paste0(outdir,"/subclass/all/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".all_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top100 hyper DHMRs position of all subclasses ##
if (!dir.exists(paste0(outdir,"/subclass/top100"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/subclass/top100"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$cluster)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(cluster) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:100)
    df = data.frame(chr=top_df[top_df$cluster == cl,"chr"],start=top_df[top_df$cluster == cl,"start"],end=top_df[top_df$cluster == cl,"end"])
    write.table(df,file=paste0(outdir,"/subclass/top100/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top100_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top500 hyper DHMRs position of all subclasses ##
if (!dir.exists(paste0(outdir,"/subclass/top500"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/subclass/top500"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$cluster)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(cluster) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:500)
    df = data.frame(chr=top_df[top_df$cluster == cl,"chr"],start=top_df[top_df$cluster == cl,"start"],end=top_df[top_df$cluster == cl,"end"])
    write.table(df,file=paste0(outdir,"/subclass/top500/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top500_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top1000 hyper DHMRs position of all subclasses ##
if (!dir.exists(paste0(outdir,"/subclass/top1000")) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/subclass/top1000")  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$cluster)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(cluster) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:1000)
    df = data.frame(chr=top_df[top_df$cluster == cl,"chr"],start=top_df[top_df$cluster == cl,"start"],end=top_df[top_df$cluster == cl,"end"])
    write.table(df,file=paste0(outdir,"/subclass/top1000/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top1000_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## all hyper DHMRs position of three classes ##
if (!dir.exists(paste0(outdir,"/three_class/all"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/all"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    df = data.frame(chr=DHMR_result.sig.sorted[DHMR_result.sig.sorted$three_class == cl,"chr"],start=DHMR_result.sig.sorted[DHMR_result.sig.sorted$three_class == cl,"start"],end=DHMR_result.sig.sorted[DHMR_result.sig.sorted$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"/three_class/all/",cl,".all_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top100 hyper DHMRs position of three classes ##
if (!dir.exists(paste0(outdir,"/three_class/top100"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/top100"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(three_class) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:100)
    df = data.frame(chr=top_df[top_df$three_class == cl,"chr"],start=top_df[top_df$three_class == cl,"start"],end=top_df[top_df$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"/three_class/top100/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top100_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top500 hyper DHMRs position of three classes ##
if (!dir.exists(paste0(outdir,"/three_class/top500"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/top500"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(three_class) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:500)
    df = data.frame(chr=top_df[top_df$three_class == cl,"chr"],start=top_df[top_df$three_class == cl,"start"],end=top_df[top_df$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"/three_class/top500/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top500_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## top1000 hyper DHMRs position of three classes ##
if (!dir.exists(paste0(outdir,"/three_class/top1000"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/top1000"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    top_df <- DHMR_result.sig.sorted %>%
        group_by(three_class) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:1000)
    df = data.frame(chr=top_df[top_df$three_class == cl,"chr"],start=top_df[top_df$three_class == cl,"start"],end=top_df[top_df$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"/three_class/top1000/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top1000_hyper_DHMRs.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

## generate mm10 gene bed ##
gene.df = read.table("../../../04.data/01.ref/mm10.genes_duplicated.bed",sep = "\t",header=T)
head(gene.df)
write.table(gene.df[,1:4],file="../../../output/03.Aging_Mouse/mm10_gene.bed",quote=F,col.names=F,row.names=F,sep="\t")

## top500 hyper DHMRs position of all subclassed of each three class ##
if (!dir.exists(paste0(outdir,"/three_class/top500"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/top500"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    DHMR_result.sig.sorted.subset = subset(DHMR_result.sig.sorted,three_class == cl)
    top_df <- DHMR_result.sig.sorted.subset %>%
        group_by(cluster) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:500)
    df = data.frame(chr=top_df[top_df$three_class == cl,"chr"],start=top_df[top_df$three_class == cl,"start"],end=top_df[top_df$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"three_class/top500/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top500_hyper_DHMRs_per_subclass.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}

head(DHMR_result.sig.sorted)

## top1000 hyper DHMRs position of all subclassed of each three class ##
if (!dir.exists(paste0(outdir,"/three_class/top1000"))) {  
            # 如果文件夹不存在，则创建它  
  dir.create(paste0(outdir,"/three_class/top1000"))  
            #cat("文件夹 'my_new_folder' 已成功创建。\n")  
}

for(cl in unique(DHMR_result.sig.sorted$three_class)){
    DHMR_result.sig.sorted.subset = subset(DHMR_result.sig.sorted,three_class == cl)
    top_df <- DHMR_result.sig.sorted.subset %>%
        group_by(cluster) %>%
        arrange(desc(diff)) %>%         # 按照diff列降序排序
        slice(1:1000)
    df = data.frame(chr=top_df[top_df$three_class == cl,"chr"],start=top_df[top_df$three_class == cl,"start"],end=top_df[top_df$three_class == cl,"end"])
    write.table(df,file=paste0(outdir,"/three_class/top1000/",str_replace_all(str_replace_all(cl," ","_"),"/","_"),".top1000_hyper_DHMRs_per_subclass.bed"),quote=F,col.names=F,row.names=F,sep="\t")
}


