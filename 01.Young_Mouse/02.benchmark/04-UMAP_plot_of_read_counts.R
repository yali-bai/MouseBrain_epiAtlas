library(Seurat)
library(stringr)
library(ggplot2)

Joint_Cabernet = readRDS("Joint_Cabernet.with_celltype.rds")
umap_loc = data.frame(Embeddings(Joint_Cabernet, reduction = 'umap'))
rownames(umap_loc) = unlist(lapply(rownames(umap_loc), function(x) strsplit(x,"@@_")[[1]][2]))

RNA_stat = read.csv("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.RNA_QC_stat.young.csv",header=T)
rownames(RNA_stat) = RNA_stat$SampleID

RNA_stat$umap_1 = umap_loc[RNA_stat$SampleID,"umap_1"]
RNA_stat$umap_2 = umap_loc[RNA_stat$SampleID,"umap_2"]

inte.col = c("darkgreen","#20A486FF", "#E8E419FF","#FDE725FF")
pdf("Raw_reads.RNA.umap.pdf",width = 5,height = 4)
ggplot(RNA_stat, aes(x = umap_1, y = umap_2, color = Total_reads)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Raw reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(1e4,3e6)) #, limits = c(1e6,3e7)
dev.off()

pdf("Mapped_reads.RNA.umap.pdf",width = 5,height = 4)
ggplot(RNA_stat, aes(x = umap_1, y = umap_2, color = reads_aligned)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Mapped reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(1e4,3e6)) #, limits = c(1e6,3e7)
dev.off()

DNA_stat = read.csv("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.DNA_QC_stat.young.csv",header=T)
paired_info = read.csv("../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",header=T)
rownames(paired_info) = paired_info$hmC_SampleID
paired_info_mC = paired_info
rownames(paired_info_mC) = paired_info_mC$mC_SampleID
merge_info = rbind(paired_info,paired_info_mC)
DNA_stat$RNA_SampleID = merge_info[DNA_stat$SampleID,"RNA_SampleID"]
DNA_stat$umap_1 = umap_loc[DNA_stat$RNA_SampleID,"umap_1"]
DNA_stat$umap_2 = umap_loc[DNA_stat$RNA_SampleID,"umap_2"]

hmC_stat = DNA_stat[DNA_stat$total_QC == 1 & DNA_stat$Library == "hmC",]
dim(hmC_stat)
mC_stat = DNA_stat[DNA_stat$total_QC == 1 & DNA_stat$Library == "mC",]
dim(mC_stat)

hmC_stat$Total_reads[hmC_stat$Total_reads > 5e7] = 5e7
hmC_stat$Aligned_Reads[hmC_stat$Aligned_Reads > 5e7] = 5e7
pdf("Raw_reads.hmC.umap.pdf",width = 5,height = 4)
ggplot(hmC_stat, aes(x = umap_1, y = umap_2, color = Total_reads)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Raw reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(8e4,5e7)) #, limits = c(1e6,3e7)
dev.off()

pdf("Mapped_reads.hmC.umap.pdf",width = 5,height = 4)
ggplot(hmC_stat, aes(x = umap_1, y = umap_2, color = Aligned_Reads)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Mapped reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(8e4,5e7)) #, limits = c(1e6,3e7)
dev.off()


pdf("Raw_reads.mC.umap.pdf",width = 5,height = 4)
ggplot(mC_stat, aes(x = umap_1, y = umap_2, color = Total_reads)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Raw reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(8e4,5e7)) #, limits = c(1e6,3e7)
dev.off()

pdf("Mapped_reads.mC.umap.pdf",width = 5,height = 4)
ggplot(mC_stat, aes(x = umap_1, y = umap_2, color = Aligned_Reads)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "Mapped reads") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
     scale_color_gradientn(colors = inte.col, limits = c(8e4,5e7)) #, limits = c(1e6,3e7)
dev.off()


RNA_temp = data.frame(type = c(rep("Raw",dim(RNA_stat)[1]),rep("Mapped",dim(RNA_stat)[1])),counts = c(RNA_stat$Total_reads,RNA_stat$reads_aligned))
RNA_temp$type = factor(RNA_temp$type,levels = c("Raw","Mapped"))
pdf("RNA_raw_mapped_reads_counts.violin_plot.pdf",width = 4,height = 3.5)
ggplot(RNA_temp, aes(x = type, y = counts,fill=type)) +
    geom_violin(linewidth=0.2,width=0.8,na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = c('Raw'='#628ebf',"Mapped"='#f1bc24'))+
    labs(title = "RNA") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),                  
            axis.text = element_text(size=13,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
    scale_y_continuous(limits=c(0,2e6),breaks=c(0,1e6,2e6),labels = c("0","1e6","2e6"))
dev.off()


hmC_temp = data.frame(type = c(rep("Raw",dim(hmC_stat)[1]),rep("Mapped",dim(hmC_stat)[1])),counts = c(hmC_stat$Total_reads,hmC_stat$Aligned_Reads))
hmC_temp$type = factor(hmC_temp$type,levels = c("Raw","Mapped"))
pdf("hmC_raw_mapped_reads_counts.violin_plot.pdf",width = 4,height = 3.5)
ggplot(hmC_temp, aes(x = type, y = counts,fill=type)) +
    geom_violin(linewidth=0.2,width=0.8,na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = c('Raw'='#628ebf',"Mapped"='#f1bc24'))+
    labs(title = "5hmC") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),                
            axis.text = element_text(size=13,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
    scale_y_continuous(limits=c(0,2e7),breaks=c(0,1e7,2e7),labels = c("0","1e7","2e7"))
dev.off()


mC_temp = data.frame(type = c(rep("Raw",dim(mC_stat)[1]),rep("Mapped",dim(mC_stat)[1])),counts = c(mC_stat$Total_reads,mC_stat$Aligned_Reads))
mC_temp$type = factor(mC_temp$type,levels = c("Raw","Mapped"))
pdf("mC_raw_mapped_reads_counts.violin_plot.pdf",width = 4,height = 3.5)
ggplot(mC_temp, aes(x = type, y = counts,fill=type)) +
    geom_violin(linewidth=0.2,width=0.8,na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = c('Raw'='#628ebf',"Mapped"='#f1bc24'))+
    labs(title = "5hmC+5mC") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),                  
            axis.text = element_text(size=13,face="bold",color = "black"),
            legend.position = "right",
            legend.key.size = unit(7, "pt"),
            legend.title = element_text(face="bold",size=8),
            legend.text = element_text(face="bold",size=7),
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())+
    scale_y_continuous(limits=c(0,2e7),breaks=c(0,1e7,2e7),labels=c("0","1e7","2e7"))
dev.off()