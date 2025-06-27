##### 01.import packages #####
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(ggpointdensity) 
library(RColorBrewer)
library(cowplot)
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/analysisdata/Methyl/public/rna/lib/R/library","/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library"))
library(MuDataSeurat)
library(RIdeogram)

##### 02.plot all chrom #####
chrom_info<-read.csv("../../04.data/01.ref/chrom_info_100k.csv")
chrom_info<-chrom_info[,c(1,2,4,3)]
mouse_chrom<-data.frame()
for(chr in unique(chrom_info$chrom)){
    chr_data<-chrom_info[chrom_info$chrom==chr,]
    df<-data.frame(Chr=chr,Start=0,End=chr_data$end[nrow(chr_data)])
    mouse_chrom<-rbind(mouse_chrom,df)
}
mouse_chrom$Chr<-gsub("chr","",mouse_chrom$Chr)
mouse_chrom<-mouse_chrom[order(as.numeric(mouse_chrom$Chr)),]


# setwd("./")

##### 03.plot DMRs #####
#plot_data = read.table("DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.count.bed",header=F)
plot_data = read.table("../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.count.bed",header=F)
colnames(plot_data) = c("Chr","Start","End","Value")
summary(plot_data$Value)
plot_data$Chr<-gsub("chr","",plot_data$Chr)
plot_data$Start <- as.numeric(plot_data$Start)
plot_data$End <- as.numeric(plot_data$End)
plot_data$Value[which(plot_data$Value > 300)] = 300
ideogram(karyotype = mouse_chrom, overlaid = plot_data,width=200,output=paste0("DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.svg"),colorset1=c("#f7f7f7", "#e34a33"))
svg2pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.svg",file="../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.pdf",width=5,height=7)

##### 04.plot DHMRs #####
#plot_data = read.table("DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.count.bed",header=F)
plot_data = read.table("../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.count.bed",header=F)
colnames(plot_data) = c("Chr","Start","End","Value")
summary(plot_data$Value)
#plot_data$Value = 1
plot_data$Chr<-gsub("chr","",plot_data$Chr)
plot_data$Start <- as.numeric(plot_data$Start)
plot_data$End <- as.numeric(plot_data$End)
plot_data$Value[which(plot_data$Value > 350)] = 350
ideogram(karyotype = mouse_chrom, overlaid = plot_data,width=200,output=paste0("DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.svg"),colorset1=c("#f7f7f7", "#e34a33"))
svg2pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.svg",file="../../output/03.Aging_Mouse/03-DMRs_DHMRs/chrom_density/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.total.pdf",width=5,height=7)



