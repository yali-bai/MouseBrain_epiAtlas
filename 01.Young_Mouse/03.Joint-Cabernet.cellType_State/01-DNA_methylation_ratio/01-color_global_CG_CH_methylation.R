library(Seurat)
library(stringr)
library(ggplot2)
library(cowplot)

Joint_Cabernet = readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/Joint_Cabernet.with_celltype.rds")
umap_loc = data.frame(Embeddings(Joint_Cabernet, reduction = 'umap'))
rownames(umap_loc) = unlist(lapply(rownames(umap_loc), function(x) strsplit(x,"@@_")[[1]][2]))

DNA_stat = read.csv("TSO-joint.DNA_QC_stat.young.add_celltype.csv",header=T)
paired_info = read.csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",header=T)
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

hmCH_df = read.table("../../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.5hmCH.global_methy.txt",header=F)
colnames(hmCH_df) = c("SampleID","mc","cov","fraction")
hmCH_df$SampleID = unlist(lapply(hmCH_df$SampleID, function(x) strsplit(x,"allc_")[[1]][2]))
hmCH_df$SampleID = str_replace_all(hmCH_df$SampleID,".mm10.dna.tsv.gz","")
hmCH_df[1:3,]                        
rownames(hmCH_df) = hmCH_df$SampleID
hmC_stat$hmCH_ratio = hmCH_df[hmC_stat$SampleID,"fraction"]
hmC_stat[1:3,]

mCH_hmCH_df = read.table("../../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.5mCH_5hmCH.global_methy.txt",header=F)
colnames(mCH_hmCH_df) = c("SampleID","mc","cov","fraction")
mCH_hmCH_df$SampleID = unlist(lapply(mCH_hmCH_df$SampleID, function(x) strsplit(x,"allc_")[[1]][2]))
mCH_hmCH_df$SampleID = str_replace_all(mCH_hmCH_df$SampleID,".mm10.dna.tsv.gz","")
rownames(mCH_hmCH_df) = mCH_hmCH_df$SampleID
mC_stat$mCH_hmCH_ratio = mCH_hmCH_df[mC_stat$SampleID,"fraction"]
mC_stat[1:3,]

pdf("global_hmCG.umap.pdf",width = 5,height = 4)
ggplot(hmC_stat, aes(x = umap_1, y = umap_2, color = dna_mCG_R)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5hmCG") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "yellow")
dev.off()

pdf("global_mCG_hmCG.umap.pdf",width = 5,height = 4)
ggplot(mC_stat, aes(x = umap_1, y = umap_2, color = dna_mCG_R)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5mCG+5hmCG") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "yellow",limits = c(0.3, 1))
dev.off()

merge_df = merge(hmC_stat,mC_stat,by="RNA_SampleID")
merge_df$mCG_ratio = merge_df$dna_mCG_R.y - merge_df$dna_mCG_R.x
merge_df$mCH_ratio = merge_df$mCH_hmCH_ratio - merge_df$hmCH_ratio
merge_df$mCH_ratio[merge_df$mCH_ratio < 0] = 0

pdf("global_mCG.umap.pdf",width = 5,height = 4)
ggplot(merge_df, aes(x = umap_1.x, y = umap_2.x, color = mCG_ratio)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5mCG") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "yellow",limits = c(0.3, 1))
dev.off()




pdf("global_mCH.umap.pdf",width = 5,height = 4)
ggplot(merge_df, aes(x = umap_1.x, y = umap_2.x, color = mCH_ratio)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5mCH") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "#75D054FF",limits = c(0, 0.03))
dev.off()

pdf("global_mCH_hmCH.umap.pdf",width = 5,height = 4)
ggplot(mC_stat, aes(x = umap_1, y = umap_2, color = mCH_hmCH_ratio)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5mCH+5hmCH") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "#75D054FF",limits = c(0, 0.03))
dev.off()

pdf("global_hmCH.umap.pdf",width = 5,height = 4)
ggplot(hmC_stat, aes(x = umap_1, y = umap_2, color = hmCH_ratio)) +
    geom_point(size=0.5) +
    labs(x = "UMAP_1", y = "UMAP_2", title = "5hmCH") +
    theme_minimal(base_family = "ArialMT") +
    theme(text = element_text(size = 10), panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_color_gradient(low = "darkblue",high = "#75D054FF",limits = c(0, 0.03))
dev.off()


##### violin plot #####
subclass_order = readRDS("../../../03.data/04.config_files/subclass_order.rds")
subclass_color = readRDS("../../../03.data/04.config_files/subclass_color.rds")
merge_df$subclass_label.y = factor(merge_df$subclass_label.y,levels = subclass_order)


p1 = ggplot(data = merge_df,aes(x = subclass_label.y, y = dna_mCG_R.x,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5hmCG\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.4), breaks=c(0,0.4),labels = function(x) sprintf("%.1f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8,face="bold",color = "black")) 

p2 = ggplot(data = merge_df,aes(x = subclass_label.y, y = mCG_ratio,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5mCG\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0.3, 1), breaks=c(0.3,1),labels = function(x) sprintf("%.1f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8,face="bold",color = "black"))

p3 = ggplot(data = merge_df,aes(x = subclass_label.y, y = dna_mCG_R.y,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5mCG+5hmCG\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0.3, 1), breaks=c(0.3,1),labels = function(x) sprintf("%.1f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,size=10,face="bold",color = "black"),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.title = element_text(size=8,face="bold",color = "black"))

pdf("Global_CG_methylation.violin_plot.pdf",width = 8,height = 5)
plot_grid(p1,p2,p3,ncol=1,rel_heights=c(1,1,2.5),align = "v",         
          axis = "l",           
          greedy = FALSE,       
          vjust = 0,           
          scale = 0.95)
dev.off()


p1 = ggplot(data = merge_df,aes(x = subclass_label.y, y = hmCH_ratio,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5hmCH\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.03), breaks=c(0, 0.03),labels = function(x) sprintf("%.2f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8,face="bold",color = "black"))

p2 = ggplot(data = merge_df,aes(x = subclass_label.y, y = mCH_ratio,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5mCH\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.03), breaks=c(0, 0.03),labels = function(x) sprintf("%.2f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8,face="bold",color = "black"))

p3 = ggplot(data = merge_df,aes(x = subclass_label.y, y = mCH_hmCH_ratio,fill = subclass_label.y)) +
    geom_violin(linewidth=0.2,width=0.8,scale="width",position = position_dodge(width = 0.4),na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = subclass_color)+
    labs(y = paste0("5mCH+5hmCH\n(Mean Frac.)"),x = "Subclass") +
    theme_bw()+
    scale_y_continuous(limits=c(0, 0.03), breaks=c(0, 0.03),labels = function(x) sprintf("%.2f", x))+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,size=10,face="bold",color = "black"),
            axis.text.y = element_text(size=10,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            axis.line=element_line(linewidth=0.5),
            plot.title = element_text(size = 15,face="bold",hjust = 0.5),
            axis.title = element_text(size=8,face="bold",color = "black"))

pdf("Global_CH_methylation.violin_plot.pdf",width = 8,height = 5)
plot_grid(p1,p2,p3,ncol=1,rel_heights=c(1,1,2.5),align = "v",          
          axis = "l",           
          greedy = FALSE,       
          vjust = 0,            
          scale = 0.95)
dev.off()
