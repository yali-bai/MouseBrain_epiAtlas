##### 01.import packages #####
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/analysisdata/Methyl/public/rna/lib/R/library","/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library"))
library(data.table)
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
# setwd("./")

##### 03.read significant DMRs DHMRs info in #####
## segment and subclass information ##
DMR_result.sig = readRDS("../../output/03.Aging_Mouse/03-DMRs_DHMRs/DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")
DHMR_result.sig=readRDS("../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")

## subclass mean 5hmCG of all segments ##
hmC_methy = fread(paste0(indir,"/aging_DMR_DHMR_5hmCG_mean_methy_level_of_subclass_diff_age.csv"),data.table=F,header=T,drop=1)
hmC_methy$type = "5hmCG"
head(hmC_methy)

## subclass mean 5mCG_5hmCG of all segments ##
mC_hmC_methy = fread(paste0(indir,"/aging_DMR_DHMR_5mCG_5hmCG_mean_methy_level_of_subclass_diff_age.csv"),data.table=F,header=T,drop=1)
mC_hmC_methy$type = "5mCG_5hmCG"

## subclass mean 5mCG of all segments ##
mC_methy = fread(paste0(indir,"/aging_DMR_DHMR_true_5mCG_mean_methy_level_of_subclass_diff_age.csv"),data.table=F,header=T,drop=1)
mC_methy$type = "5mCG"

## merge ##
merge.df = rbind(hmC_methy,mC_hmC_methy,mC_methy)
merge.df$subclass = unlist(lapply(as.character(merge.df$unique), function(x) strsplit(x,'_')[[1]][1]))
merge.df$age = unlist(lapply(as.character(merge.df$unique), function(x) strsplit(x,'_')[[1]][2]))
unique(merge.df$subclass)
unique(merge.df$age)

## deformation for calculating mean difference ##
merge.df$unique = paste0(merge.df$subclas,";",merge.df$segment,";",merge.df$type)
merge.df.dcast = dcast(merge.df,unique~age)
head(merge.df.dcast)
merge.df.dcast$diff = merge.df.dcast$old - merge.df.dcast$young

merge.df.dcast$subclass = unlist(lapply(as.character(merge.df.dcast$unique), function(x) strsplit(x,';')[[1]][1]))
merge.df.dcast$segment = unlist(lapply(as.character(merge.df.dcast$unique), function(x) strsplit(x,';')[[1]][2]))
merge.df.dcast$type = unlist(lapply(as.character(merge.df.dcast$unique), function(x) strsplit(x,';')[[1]][3]))
head(merge.df.dcast)                                       

## save result for plotting heatmap ##
saveRDS(merge.df.dcast,file=paste0(outdir,"/merge.df.dcast.rds"))


## subclass order ##
subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")

## hyper DMRs ## 
## generate data for boxplot ##
DMR_plot.df = data.frame()
for (cl in unique(DMR_result.sig$cluster)){
    seg.v = DMR_result.sig[DMR_result.sig$cluster == cl & DMR_result.sig$diff > 0.05,'chrom']
    idx = intersect(which(merge.df.dcast$subclass == cl),which(merge.df.dcast$segment %in% seg.v))
    DMR_plot.df = rbind(DMR_plot.df,merge.df.dcast[idx,])  
}
DMR_plot.df$subclass = factor(DMR_plot.df$subclass,levels=subclass_order)
p1<-ggplot(DMR_plot.df,aes(x=subclass,y=diff,fill=type))+
    geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
    ylab("mean diff (old - young)")+
    xlab("subclass")+
    scale_fill_manual(values = c('5hmCG'='#3498db','5mCG'='#e74c3c','5mCG_5hmCG'="#ffdc75"))+
    theme_bw()+
    coord_cartesian(ylim = c(-0.5, 0.5))+
    geom_hline(yintercept=0, linetype='dashed', color='red')+
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
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/significant_DMRs_DHMRs_boxplot/aging_hyper_DMR_DNA_methyl_level_boxplot.pdf",width=12,height=3)   
print(p1)
dev.off()

## hypo DMRs ## 
## generate data for boxplot ##
DMR_plot.df = data.frame()
for (cl in unique(DMR_result.sig$cluster)){
    seg.v = DMR_result.sig[DMR_result.sig$cluster == cl & DMR_result.sig$diff < -0.05,'chrom']
    idx = intersect(which(merge.df.dcast$subclass == cl),which(merge.df.dcast$segment %in% seg.v))
    DMR_plot.df = rbind(DMR_plot.df,merge.df.dcast[idx,])  
}
DMR_plot.df$subclass = factor(DMR_plot.df$subclass,levels=subclass_order)
p1<-ggplot(DMR_plot.df,aes(x=subclass,y=diff,fill=type))+
    geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
    ylab("mean diff (old - young)")+
    xlab("subclass")+
    scale_fill_manual(values = c('5hmCG'='#3498db','5mCG'='#e74c3c','5mCG_5hmCG'="#ffdc75"))+
    theme_bw()+
    coord_cartesian(ylim = c(-0.75, 0.5))+
    geom_hline(yintercept=0, linetype='dashed', color='red')+
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
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/significant_DMRs_DHMRs_boxplot/aging_hypo_DMR_DNA_methyl_level_boxplot.pdf",width=10,height=3)   
print(p1)
dev.off()

## hyper DHMRs ## 
## generate data for boxplot ##
DHMR_plot.df = data.frame()
for (cl in unique(DHMR_result.sig$cluster)){
    seg.v = DHMR_result.sig[DHMR_result.sig$cluster == cl & DHMR_result.sig$diff > 0.05,'chrom']
    idx = intersect(which(merge.df.dcast$subclass == cl),which(merge.df.dcast$segment %in% seg.v))
    DHMR_plot.df = rbind(DHMR_plot.df,merge.df.dcast[idx,])  
}
DHMR_plot.df$subclass = factor(DHMR_plot.df$subclass,levels=subclass_order)
p1<-ggplot(DHMR_plot.df,aes(x=subclass,y=diff,fill=type))+
    geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
    ylab("mean diff (old - young)")+
    xlab("subclass")+
    scale_fill_manual(values = c('5hmCG'='#3498db','5mCG'='#e74c3c','5mCG_5hmCG'="#ffdc75"))+
    theme_bw()+
    coord_cartesian(ylim = c(-0.5, 0.5))+
    geom_hline(yintercept=0, linetype='dashed', color='red')+
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
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/significant_DMRs_DHMRs_boxplot/aging_hyper_DHMR_DNA_methyl_level_boxplot.pdf",width=12,height=3)   
print(p1)
dev.off()

## hypO DHMRs ## 
## generate data for boxplot ##
DHMR_plot.df = data.frame()
for (cl in unique(DHMR_result.sig$cluster)){
    seg.v = DHMR_result.sig[DHMR_result.sig$cluster == cl & DHMR_result.sig$diff < -0.05,'chrom']
    idx = intersect(which(merge.df.dcast$subclass == cl),which(merge.df.dcast$segment %in% seg.v))
    DHMR_plot.df = rbind(DHMR_plot.df,merge.df.dcast[idx,])
}
DHMR_plot.df$subclass = factor(DHMR_plot.df$subclass,levels=subclass_order)
p1<-ggplot(DHMR_plot.df,aes(x=subclass,y=diff,fill=type))+
    geom_boxplot(na.rm = TRUE,outlier.shape = NA)+
    ylab("mean diff (old - young)")+
    xlab("subclass")+
    scale_fill_manual(values = c('5hmCG'='#3498db','5mCG'='#e74c3c','5mCG_5hmCG'="#ffdc75"))+
    theme_bw()+
    coord_cartesian(ylim = c(-0.5, 0.5))+
    geom_hline(yintercept=0, linetype='dashed', color='red')+
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
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))
pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/significant_DMRs_DHMRs_boxplot/aging_hypo_DHMR_DNA_methyl_level_boxplot.pdf",width=12,height=3)
print(p1)
dev.off()
