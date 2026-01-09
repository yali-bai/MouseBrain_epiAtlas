##### import packages #####
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(stringr)

outdir = ""

DNA_stat = read.csv("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.DNA_QC_stat.young.csv",header=T)
head(DNA_stat)

hmC_stat = DNA_stat[DNA_stat$total_QC == 1 & DNA_stat$Library == "hmC",]
dim(hmC_stat)
mC_stat = DNA_stat[DNA_stat$total_QC == 1 & DNA_stat$Library == "mC",]
dim(mC_stat)

DNA_stat[1,"dna_mCG_R"]

pdf(paste0(outdir,"/","hmC.dna_mCG_R.histogram.pdf"),width=5,height=1.8)
hmC_stat$dna_mCG_R = as.numeric(hmC_stat$dna_mCG_R)*100
ggplot(hmC_stat, aes(x=dna_mCG_R))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("hmCG(%)")+
    scale_y_continuous(limits=c(0, 4000), breaks=c(0,2000,4000))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()
pdf(paste0(outdir,"/","mC.dna_mCG_R.histogram.pdf"),width=5,height=1.8)
mC_stat$dna_mCG_R = as.numeric(mC_stat$dna_mCG_R)*100
ggplot(mC_stat, aes(x=dna_mCG_R))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("mCG+hmCG(%)")+
    scale_y_continuous(limits=c(0, 4000), breaks=c(0,2000,4000))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()

pdf(paste0(outdir,"/","hmC.fullpuc19_mCG_R_after.histogram.pdf"),width=5,height=1.8)
hmC_stat$fullpuc19_mCG_R_after = as.numeric(hmC_stat$fullpuc19_mCG_R_after)*100
ggplot(hmC_stat, aes(x=fullpuc19_mCG_R_after))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("pUC19(%, Control)")+
    scale_y_continuous(limits=c(0, 5500), breaks=c(0,2500,5000))+
    scale_x_continuous(limits=c(0, 3), breaks=c(0,1,2,3))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()
pdf(paste0(outdir,"/","mC.fullpuc19_mCG_R_after.histogram.pdf"),width=5,height=1.8)
mC_stat$fullpuc19_mCG_R_after = as.numeric(mC_stat$fullpuc19_mCG_R_after)*100
ggplot(mC_stat, aes(x=fullpuc19_mCG_R_after))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("pUC19(%, Control)")+
    scale_y_continuous(limits=c(0, 5500), breaks=c(0,2500,5000))+
    scale_x_continuous(limits=c(96, 100), breaks=c(96,97,98,99,100))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()

pdf(paste0(outdir,"/","hmC.lambda_mCG_R_after.histogram.pdf"),width=5,height=1.8)
ggplot(hmC_stat, aes(x=lambda_mCG_R_after))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("Lambda DNA(%, Control)")+
    scale_y_continuous(limits=c(0, 10000), breaks=c(0,5000,10000))+
    scale_x_continuous(limits=c(0, 2), breaks=c(0,1,2))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()
pdf(paste0(outdir,"/","mC.lambda_mCG_R_after.histogram.pdf"),width=5,height=1.8)
ggplot(mC_stat, aes(x=lambda_mCG_R_after))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("Lambda DNA(%, Control)")+
    scale_y_continuous(limits=c(0, 10000), breaks=c(0,5000,10000))+
    scale_x_continuous(limits=c(0, 2), breaks=c(0,1,2))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()

pdf(paste0(outdir,"/","hmC.cla1_mCG_R_pre.histogram.pdf"),width=5,height=1.8)
ggplot(hmC_stat, aes(x=cla1_mCG_R_pre))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("Cla1(%, Control)")+
    scale_y_continuous(limits=c(0, 5000), breaks=c(0,2500,5000))+
    scale_x_continuous(limits=c(98, 100), breaks=c(98,99,100))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()
pdf(paste0(outdir,"/","mC.cla1_mCG_R_pre.histogram.pdf"),width=5,height=1.8)
ggplot(mC_stat, aes(x=cla1_mCG_R_pre))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("Cla1(%, Control)")+
    scale_y_continuous(limits=c(0, 5000), breaks=c(0,2500,5000))+
    scale_x_continuous(limits=c(97, 100), breaks=c(97,98,99,100))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()

hmCH_df = read.table("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.5hmCH.global_methy.txt",header=F)
hmCH_df[1:3,]

colnames(hmCH_df) = c("SampleID","mc","cov","fraction")
hmCH_df$SampleID = unlist(lapply(hmCH_df$SampleID, function(x) strsplit(x,"allc_")[[1]][2]))
hmCH_df$SampleID = str_replace_all(hmCH_df$SampleID,".mm10.dna.tsv.gz","")
hmCH_df[1:3,]                        

rownames(hmCH_df) = hmCH_df$SampleID
hmC_stat$hmCH_ratio = hmCH_df[hmC_stat$SampleID,"fraction"]
hmC_stat[1:3,]

mCH_hmCH_df = read.table("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.5mCH_5hmCH.global_methy.txt",header=F)
colnames(mCH_hmCH_df) = c("SampleID","mc","cov","fraction")
mCH_hmCH_df$SampleID = unlist(lapply(mCH_hmCH_df$SampleID, function(x) strsplit(x,"allc_")[[1]][2]))
mCH_hmCH_df$SampleID = str_replace_all(mCH_hmCH_df$SampleID,".mm10.dna.tsv.gz","")
rownames(mCH_hmCH_df) = mCH_hmCH_df$SampleID
mC_stat$mCH_hmCH_ratio = mCH_hmCH_df[mC_stat$SampleID,"fraction"]
mC_stat[1:3,]

pdf(paste0(outdir,"/","hmCH.histogram.pdf"),width=5,height=1.8)
ggplot(hmC_stat, aes(x=hmCH_ratio))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("hmCH(%)")+
    scale_y_continuous(limits=c(0, 10000), breaks=c(0,5000,10000))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()


pdf(paste0(outdir,"/","mCH_hmCH.histogram.pdf"),width=5,height=1.8)
ggplot(mC_stat, aes(x=mCH_hmCH_ratio))+
    geom_histogram(linewidth=0.2,color = "black",fill="#15869d", position = "dodge",bins = 60,width = 0.5)+
    ggtitle("mCH+hmCH(%)")+
    scale_y_continuous(limits=c(0, 10000), breaks=c(0,5000,10000))+
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black",face="bold", size=15),
        axis.title = element_blank(),
        plot.title = element_text(size = 18,face = "bold",color = "black")
    )
dev.off()


