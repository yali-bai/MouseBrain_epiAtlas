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
library(ggunchained)
library(reshape2)
library(ggpubr)
library(cowplot)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.change working path #####
setwd('./')
##### 03. data prepare #####
## load integrated data and metainfo ##
load(paste0(indir,"/seurat_obj.Joint_Cabernet_integrated.RData"))
metainfo = our.seuratobj@meta.data

## get unique id for pairing ##
metainfo$sample = merged.seuratobj.sct@meta.data[rownames(metainfo),"sample"]
metainfo$uniq_id = paste(unlist(lapply(metainfo$sample, function(x) strsplit(x,"_")[[1]][1])),unlist(lapply(metainfo$sample, function(x) strsplit(x,"_")[[1]][2])),unlist(lapply(metainfo$sample, function(x) strsplit(x,"_")[[1]][4])),unlist(lapply(metainfo$sample, function(x) strsplit(x,"_")[[1]][6])),unlist(lapply(metainfo$sample, function(x) strsplit(x,"_")[[1]][8])),sep="_")

young_match.df = read.csv("../../04.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",header=T)
young_match.subset = subset(young_match.df,RNA_QC==1)
metainfo$uniq_id[match(young_match.subset$RNA,metainfo$sample)] = young_match.subset$Unique_ID_match

## QC information ##
hmC_QC.df = fread("../../04.data/02.metainfo/03.Aging_Mouse/TSO-joint.hmC_QC_stat.aged.csv",header=T,,data.table=F)
mC_QC.df = fread("../../04.data/02.metainfo/03.Aging_Mouse/TSO-joint.mC_QC_stat.aged.csv",header=T,,data.table=F)

## DNA global level ##
## CH ##
## hmCH ##
global_methy = fread("../../04.data/05.intermediate_files/02.DNA/02.Aging_Mouse/5hmCH.all_cells.global_methy.txt",header=F,data.table=F)
global_methy$uniq_id = paste(unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][1])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][2])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][4])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][6])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][8])),sep="_")
global_methy$uniq_id[match(intersect(young_match.df$hmC,global_methy$V1),global_methy$V1)] = young_match.df$Unique_ID_match[match(intersect(young_match.df$hmC,global_methy$V1),young_match.df$hmC)]
global_methy$QC = NA
global_methy$QC[match(intersect(global_methy$V1,hmC_QC.df$SampleID),global_methy$V1)] = hmC_QC.df[match(intersect(global_methy$V1,hmC_QC.df$SampleID),hmC_QC.df$SampleID),"QC"]
global_methy$V2[which(global_methy$QC == 0)] = NA
rownames(global_methy) = global_methy$uniq_id
global_methy_subset = global_methy[metainfo$uniq_id,]
metainfo$hmCH = global_methy_subset$V2
metainfo$age = factor(metainfo$age,levels=c("young","old"))

## mCH ##
global_methy = fread("../../04.data/05.intermediate_files/02.DNA/02.Aging_Mouse/5mCH_5hmCH.all_cells.global_methy.txt",header=F,data.table=F)
global_methy$uniq_id = paste(unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][1])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][2])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][4])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][6])),unlist(lapply(global_methy$V1, function(x) strsplit(x,"_")[[1]][8])),sep="_")
global_methy$uniq_id[match(intersect(young_match.df$mC,global_methy$V1),global_methy$V1)] = young_match.df$Unique_ID_match[match(intersect(young_match.df$mC,global_methy$V1),young_match.df$mC)]
global_methy$QC = NA
global_methy$QC[match(intersect(global_methy$V1,mC_QC.df$SampleID),global_methy$V1)] = mC_QC.df[match(intersect(global_methy$V1,mC_QC.df$SampleID),mC_QC.df$SampleID),"QC"]
global_methy$V2[which(global_methy$QC == 0)] = NA
rownames(global_methy) = global_methy$uniq_id
global_methy_subset = global_methy[metainfo$uniq_id,]
metainfo$mCH_hmCH = global_methy_subset$V2
metainfo$age = factor(metainfo$age,levels=c("young","old"))

## mCH ##
metainfo$mCH_hmCH = as.numeric(metainfo$mCH_hmCH)
metainfo$hmCH = as.numeric(metainfo$hmCH)

metainfo$mCH = NA
metainfo$mCH = metainfo$mCH_hmCH - metainfo$hmCH

## CG in stat file ##
## hmCG ##
metainfo$stat_hmCG = NA
young.DNA.stat.df = fread("../../04.data/02.metainfo/01.Young_Mouse/TSO-joint.hmC_QC_stat.young.csv",header=T,,data.table=F)
young.DNA.stat.hmC = subset(young.DNA.stat.df,Library == "hmC")
rownames(young.DNA.stat.hmC) = young.DNA.stat.hmC$Unique_ID
young.DNA.stat.hmC[match(intersect(young.DNA.stat.hmC$SampleID,hmC_QC.df$SampleID[which(hmC_QC.df$QC == 0)]),young.DNA.stat.hmC$SampleID),"dna_mCG_R"] = NA
metainfo$stat_hmCG[match(intersect(metainfo$uniq_id,young.DNA.stat.hmC$Unique_ID),metainfo$uniq_id)] = young.DNA.stat.hmC$dna_mCG_R[match(intersect(metainfo$uniq_id,young.DNA.stat.hmC$Unique_ID),young.DNA.stat.hmC$Unique_ID)]
old.DNA.stat.df = read.csv("../../04.data/02.metainfo/03.Aging_Mouse/TSO-joint.hmC_QC_stat.aged.csv",header=T)
old.DNA.stat.df$Unique_ID = paste(unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][1])),unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][2])),unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][4])),unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][6])),unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][8])),sep="_")
old.DNA.stat.df$Library = str_replace(unlist(lapply(old.DNA.stat.df$SampleID, function(x) strsplit(x,"_")[[1]][7])),"joint5","")
old.DNA.stat.hmC = subset(old.DNA.stat.df,Library == "hmC")
rownames(old.DNA.stat.hmC) = old.DNA.stat.hmC$Unique_ID
old.DNA.stat.hmC[match(intersect(old.DNA.stat.hmC$SampleID,hmC_QC.df$SampleID[which(hmC_QC.df$QC == 0)]),old.DNA.stat.hmC$SampleID),"dna_mCG_R"] = NA
metainfo$stat_hmCG[match(intersect(metainfo$uniq_id,old.DNA.stat.hmC$Unique_ID),metainfo$uniq_id)] = old.DNA.stat.hmC$dna_mCG_R[match(intersect(metainfo$uniq_id,old.DNA.stat.hmC$Unique_ID),old.DNA.stat.hmC$Unique_ID)]

## mCG_hmCG ##
metainfo$stat_mCG_hmCG = NA
young.DNA.stat.mC = subset(young.DNA.stat.df,Library == "mC")
rownames(young.DNA.stat.mC) = young.DNA.stat.mC$Unique_ID
metainfo$stat_mCG_hmCG[match(intersect(metainfo$uniq_id,young.DNA.stat.mC$Unique_ID),metainfo$uniq_id)] = young.DNA.stat.mC$dna_mCG_R[match(intersect(metainfo$uniq_id,young.DNA.stat.mC$Unique_ID),young.DNA.stat.mC$Unique_ID)]
young.DNA.stat.hmC[match(intersect(young.DNA.stat.hmC$SampleID,mC_QC.df$SampleID[which(mC_QC.df$QC == 0)]),young.DNA.stat.hmC$SampleID),"dna_mCG_R"] = NA
old.DNA.stat.mC = subset(old.DNA.stat.df,Library == "mC")
rownames(old.DNA.stat.mC) = old.DNA.stat.mC$Unique_ID
old.DNA.stat.hmC[match(intersect(old.DNA.stat.hmC$SampleID,mC_QC.df$SampleID[which(mC_QC.df$QC == 0)]),old.DNA.stat.hmC$SampleID),"dna_mCG_R"] = NA
metainfo$stat_mCG_hmCG[match(intersect(metainfo$uniq_id,old.DNA.stat.mC$Unique_ID),metainfo$uniq_id)] = old.DNA.stat.mC$dna_mCG_R[match(intersect(metainfo$uniq_id,old.DNA.stat.mC$Unique_ID),old.DNA.stat.mC$Unique_ID)]

## true mCG ##
metainfo$stat_mCG_hmCG = as.numeric(metainfo$stat_mCG_hmCG)
metainfo$stat_hmCG = as.numeric(metainfo$stat_hmCG)

metainfo$stat_mCG = NA
metainfo$stat_mCG = metainfo$stat_mCG_hmCG - metainfo$stat_hmCG

## save result ##
#saveRDS(metainfo,file="metainfo.250115.rds")
saveRDS(metainfo,file="../../output/03.Aging_Mouse/metainfo.250115.rds")

## subclass order ##
subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")
metainfo$lt_twice_subclass = factor(metainfo$lt_twice_subclass,levels=subclass_order)

##### 04. plot #####
## every plot has specific data range limit ##
## subclass without significance test ##
## hmCG ##
p1<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.39),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.44), breaks=c(0,0.2,0.4),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCG ##
p2<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_mCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.78),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0.4, 0.84),breaks=c(0.4,0.6,0.8),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCG_hmCG ##
p3<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_mCG_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG_hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.9),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0.7, 0.94), breaks=c(0.7,0.8,0.9),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## hmCH ##
p4<-ggplot(metainfo,aes(x=lt_twice_subclass,y=hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.011),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0.004, 0.014), breaks=c(0.006,0.009,0.012),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCH ##
p5<-ggplot(metainfo,aes(x=lt_twice_subclass,y=mCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.038),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.044), breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCH_hmCH ##
p6<-ggplot(metainfo,aes(x=lt_twice_subclass,y=mCH_hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    scale_x_discrete("subclass")+
    ylab("mCH_hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.048),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.054), breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,color = "black",size=7),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.line = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

pdf("../../output/03.Aging_Mouse/01-global_DNA_methyl_old_vs_young_violinplot/subclass.DNA_Global_methylation.split_violin_plot.with_significance.scale_width.wilcox_test.250115.pdf",width=8,height=6)
plot_grid(p1,p2,p3,p4,p5,p6,ncol=1,rel_heights=c(1,1,1,1,1,2.5))
dev.off()

## subclass without significance test ##
## hmCG ##
metainfo$lt_twice_subclass = factor(metainfo$lt_twice_subclass,levels=subclass_order)
p1<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.44), breaks=c(0,0.2,0.4),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCG ##
p2<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_mCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0.4, 0.84),breaks=c(0.4,0.6,0.8),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

p3<-ggplot(metainfo,aes(x=lt_twice_subclass,y=stat_mCG_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG_hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0.7, 0.94), breaks=c(0.7,0.8,0.9),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## hmCH ##
p4<-ggplot(metainfo,aes(x=lt_twice_subclass,y=hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0.004, 0.014), breaks=c(0.006,0.009,0.012),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCH ##
p5<-ggplot(metainfo,aes(x=lt_twice_subclass,y=mCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.044), breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

## mCH_hmCH ##
p6<-ggplot(metainfo,aes(x=lt_twice_subclass,y=mCH_hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    scale_x_discrete("subclass")+
    ylab("mCH_hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.054), breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,color = "black",size=7),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.line = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))

pdf("../../output/03.Aging_Mouse/01-global_DNA_methyl_old_vs_young_violinplot/subclass.DNA_Global_methylation.split_violin_plot.without_significance.scale_width.wilcox_test.250115.pdf",width=8,height=6)
plot_grid(p1,p2,p3,p4,p5,p6,ncol=1,rel_heights=c(1,1,1,1,1,2.5))
dev.off()


## major class with significance test ##
## major class order ##
order.new = readRDS("../../04.data/04.config_files/order.majorclass.rds")
metainfo$lt_twice_class = factor(metainfo$lt_twice_class,levels=order.new)
## hmCG ##
p1<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.39),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.44),breaks=c(0,0.2,0.4),labels = function(x) sprintf("%.3f", x))

## mCG ##
p2<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_mCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.78),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.4, 0.84),breaks=c(0.4,0.6,0.8),labels = function(x) sprintf("%.3f", x))                       

## mCG_hmCG ##
p3<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_mCG_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG_hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.9),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.7, 0.94),breaks=c(0.7,0.8,0.9),labels = function(x) sprintf("%.3f", x))                       
                       
## hmCH ##
p4<-ggplot(metainfo,aes(x=lt_twice_class,y=hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.012),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.004, 0.014),breaks=c(0.004,0.008,0.012),labels = function(x) sprintf("%.3f", x))                       
                       
## mCH ##
p5<-ggplot(metainfo,aes(x=lt_twice_class,y=mCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.044),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.049),breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))                       
                                        
## mCH_hmCH ##
p6<-ggplot(metainfo,aes(x=lt_twice_class,y=mCH_hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    scale_x_discrete("subclass")+
    ylab("mCH_hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    stat_compare_means(aes(group=age),
                    method = "wilcox.test",
                    paired = F,
                    symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                       symbols=c("***","**","*","ns")),
                    label.y = c(0.055),
                    label = "p.signif",
                    size=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,color = "black",size=7),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.line = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.06),breaks=c(0.02,0.04,0.06),labels = function(x) sprintf("%.3f", x))
pdf("../../output/03.Aging_Mouse/01-global_DNA_methyl_old_vs_young_violinplot/majorclass.DNA_Global_methylation.split_violin_plot.with_significance.scale_width.wilcox_test.20250115.pdf",width=4.5,height=6)
plot_grid(p1,p2,p3,p4,p5,p6,ncol=1,rel_heights=c(1,1,1,1,1,2))
dev.off()

## major class with significance test ##
## hmCG ##
p1<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.44),breaks=c(0,0.2,0.4),labels = function(x) sprintf("%.3f", x))

## mCG ##
p2<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_mCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.4, 0.84),breaks=c(0.4,0.6,0.8),labels = function(x) sprintf("%.3f", x))                       

## mCG_hmCG ##
p3<-ggplot(metainfo,aes(x=lt_twice_class,y=stat_mCG_hmCG,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCG_hmCG")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.7, 0.94),breaks=c(0.7,0.8,0.9),labels = function(x) sprintf("%.3f", x))                       

## hmCH ##
p4<-ggplot(metainfo,aes(x=lt_twice_class,y=hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0.004, 0.014),breaks=c(0.004,0.008,0.012),labels = function(x) sprintf("%.3f", x))                       

## mCH ##
p5<-ggplot(metainfo,aes(x=lt_twice_class,y=mCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    ylab("mCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.049),breaks=c(0,0.02,0.04),labels = function(x) sprintf("%.3f", x))                       

## mCH_hmCH ##
p6<-ggplot(metainfo,aes(x=lt_twice_class,y=mCH_hmCH,fill=age))+
    geom_split_violin(linewidth=0.2,scale = "width")+
    scale_x_discrete("subclass")+
    ylab("mCH_hmCH")+
    scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,color = "black",size=7),
            axis.text.y = element_text(size=7,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            axis.text = element_text(face="bold", size=10),
            text = element_text(face="bold",size = 10),
            axis.line = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=10))+
    scale_y_continuous(limits=c(0, 0.06),breaks=c(0.02,0.04,0.06),labels = function(x) sprintf("%.3f", x))
pdf("../../output/03.Aging_Mouse/01-global_DNA_methyl_old_vs_young_violinplot/majorclass.DNA_Global_methylation.split_violin_plot.without_significance.scale_width.wilcox_test.20250115.pdf",width=4.5,height=6)
plot_grid(p1,p2,p3,p4,p5,p6,ncol=1,rel_heights=c(1,1,1,1,1,2))
dev.off()

##### 05. calculate subclass logFC #####
## hmCG ##
hmCG_result <- metainfo %>%
  group_by(age, lt_twice_subclass) %>%
  summarise(mean_hmCG = mean(stat_hmCG, na.rm = TRUE))
head(hmCG_result)

## mCG ##
mCG_result <- metainfo %>%
  group_by(age, lt_twice_subclass) %>%
  summarise(mean_mCG = mean(stat_mCG, na.rm = TRUE))

## mCG_hmCG ##
mCG_hmCG_result <- metainfo %>%
  group_by(age, lt_twice_subclass) %>%
  summarise(mean_mCG_hmCG = mean(stat_mCG_hmCG, na.rm = TRUE))

## hmCG ##
hmCG.dcast = dcast(hmCG_result,lt_twice_subclass~age)
hmCG.dcast$hmCG_fold_change = hmCG.dcast$old/hmCG.dcast$young
head(hmCG.dcast)

## mCG ##
mCG.dcast = dcast(mCG_result,lt_twice_subclass~age)
mCG.dcast$mCG_fold_change = mCG.dcast$old/mCG.dcast$young

## mCG_hmCG ##
mCG_hmCG.dcast = dcast(mCG_hmCG_result,lt_twice_subclass~age)
mCG_hmCG.dcast$mCG_hmCG_fold_change = mCG_hmCG.dcast$old/mCG_hmCG.dcast$young

## combine data ##
merge.df = merge(hmCG.dcast,mCG.dcast,by="lt_twice_subclass")
merge.dataf = merge(merge.df,mCG_hmCG.dcast,by="lt_twice_subclass")
head(merge.dataf)
merge.reorder = merge.dataf[,c("lt_twice_subclass","hmCG_fold_change","mCG_fold_change","mCG_hmCG_fold_change")]
head(merge.reorder)

## sort by subclass order of plot ##
subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")
merge.reorder$lt_twice_subclass = factor(merge.reorder$lt_twice_subclass,levels=subclass_order)
merge.reorder = merge.reorder[order(merge.reorder$lt_twice_subclass),]
write.csv(merge.reorder,file="../../output/03.Aging_Mouse/global_DNA_old_vs_young_folg_change.csv",quote=F,row.names=F,col.names=T)


