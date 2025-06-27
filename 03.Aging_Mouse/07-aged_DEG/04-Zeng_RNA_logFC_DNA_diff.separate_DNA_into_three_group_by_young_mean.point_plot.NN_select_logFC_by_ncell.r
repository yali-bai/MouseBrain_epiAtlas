##### 01. import packages #####
library(data.table)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(ggpointdensity) 
library(RColorBrewer)
library(cowplot)
library(tidyr)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 03. Zeng DEG info #####
RNA_data = read.table(paste0(indir,"/Zeng_DEG.subclass.txt",header=T,sep="\t"))
head(RNA_data)

Zeng_NN_DEG = read.table(paste0(indir,"/Zeng_supercluster.NN.DEG.txt",header=T,sep="\t"))
head(Zeng_NN_DEG)

## change subclass names ##
temp = intersect(Zeng_NN_DEG[Zeng_NN_DEG$supertype_label == "Astro-TE NN_1","gene"],Zeng_NN_DEG[Zeng_NN_DEG$supertype_label == "Astro-TE NN_3","gene"])
length(intersect(temp,Zeng_NN_DEG[Zeng_NN_DEG$supertype_label == "Astro-TE NN_5","gene"]))

Zeng_NN_DEG$subclass_label = Zeng_NN_DEG$supertype_label
Zeng_NN_DEG$subclass_label[which(!is.na(str_match(Zeng_NN_DEG$supertype_label,"Astro-TE NN")))] = "Astro-TE NN"
Zeng_NN_DEG$subclass_label[which(Zeng_NN_DEG$supertype_label == "COP NN_1")] = 'OPC NN'
Zeng_NN_DEG$subclass_label[which(Zeng_NN_DEG$supertype_label == 'OPC NN_1')] = 'OPC NN'
Zeng_NN_DEG$subclass_label[which(Zeng_NN_DEG$supertype_label == 'MFOL NN_3')] = 'Oligo NN'
Zeng_NN_DEG$subclass_label[which(Zeng_NN_DEG$supertype_label == 'MOL NN_4')] = 'Oligo NN'
Zeng_NN_DEG$subclass_label[which(Zeng_NN_DEG$supertype_label == 'Microglia NN_1')] = 'Microglia NN'
unique(Zeng_NN_DEG$subclass_label)

RNA_data$supertype_label = NA
combined_RNA_data = rbind(RNA_data,Zeng_NN_DEG)
dim(RNA_data)
dim(combined_RNA_data )
head(combined_RNA_data )

## logFC of gene appeared in more than one supercluster is decided by ncell of supercluster ## 
supercluster_metainfo = read.table(paste0(indir,"/supercluster_metainfo.txt",header=T,sep="\t"))
head(supercluster_metainfo)

combined_RNA_data$ncell = supercluster_metainfo[match(combined_RNA_data$supertype_label,supercluster_metainfo$supertype_label),"ncell"]
tail(combined_RNA_data)

combined_RNA_data_remove_gene = combined_RNA_data %>%
    group_by(subclass_label,gene) %>%
    arrange(desc(ncell)) %>%
    slice_head(n = 1)


gene_metainfo = read.table("../../04.data/01.ref/mm10.genes_duplicated.bed",header=T)
head(gene_metainfo)
rownames(gene_metainfo)=gene_metainfo$gene_id

combined_RNA_data_remove_gene$gene_id = gene_metainfo[match(combined_RNA_data_remove_gene$gene,gene_metainfo$gene_name),"gene_id"]
head(combined_RNA_data_remove_gene)

DNA_data<-fread(paste0(indir,"/5hmC_genebody.CG_old_young_diff_result.csv"),data.table=F)
head(DNA_data)
DNA_data$gene = gene_metainfo[DNA_data$chrom,"gene_name"]
colnames(DNA_data)[colnames(DNA_data) == "cluster"]="subclass_label"
DNA_subclass<-DNA_data[grep("subclass",DNA_data$level),]
DNA_subclass<-DNA_subclass[!is.infinite(DNA_subclass$logFC)&!is.na(DNA_subclass$logFC),]
DNA_subclass$diff = DNA_subclass$old_mean - DNA_subclass$young_mean
head(DNA_subclass)

plot_data<-merge(DNA_subclass,combined_RNA_data_remove_gene,by=c("gene","subclass_label"))

freq.df = as.data.frame(table(plot_data$subclass_label))
freq.df
freq.df = freq.df[freq.df$Freq > 15,] 
freq.df

subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")

freq.df$Var1 = factor(freq.df$Var1,levels = subclass_order)
freq.df = freq.df[order(freq.df$Var1),]

head(plot_data)
plot_data_top3 = plot_data %>%
    group_by(subclass_label) %>%
    arrange(desc(age_effect_size)) %>%
    slice_head(n = 3)
plot_data_top3$yend = NA
plot_data_top3$xend = NA
for(cl in unique(plot_data_top3$subclass_label)){
    idx = which(plot_data_top3$subclass_label == cl)
    xstart = -0.03
    ystart = 0.5
    text_start = 1
    for(i in order(plot_data_top3$diff[idx])){
        plot_data_top3$yend[idx[i]] = plot_data_top3$age_effect_size[idx[i]] + ystart
        plot_data_top3$xend[idx[i]] = plot_data_top3$diff[idx[i]] + xstart
        xstart = xstart + 0.03
        ystart = ystart + 0.5
    }
}
plot_data_top3$text = plot_data_top3$yend + 0.5

plot_data_tail3 = plot_data %>%
    group_by(subclass_label) %>%
    arrange(desc(age_effect_size)) %>%
    slice_tail(n = 3)
plot_data_tail3$yend = NA
plot_data_tail3$xend = NA
for(cl in unique(plot_data_tail3$subclass_label)){
    idx = which(plot_data_tail3$subclass_label == cl)
    xstart = -0.03
    ystart = -0.5
    for(i in order(plot_data_top3$diff[idx])){
        plot_data_tail3$yend[idx[i]] = plot_data_tail3$age_effect_size[idx[i]] + ystart
        plot_data_tail3$xend[idx[i]] = plot_data_tail3$diff[idx[i]] + xstart
        xstart = xstart + 0.03
        ystart = ystart - 1
    }
}
plot_data_tail3$text = plot_data_tail3$yend - 0.5
combined_plot_data = rbind(plot_data_top3,plot_data_tail3)
head(combined_plot_data)

plot_data$diff = as.numeric(plot_data$diff)
correlation_df = data.frame()
for(sc in freq.df$Var1){
    temp = plot_data[plot_data$subclass == sc & !is.na(plot_data$age_effect_size) & !is.infinite(plot_data$age_effect_size),]
    test = cor.test(temp$diff, temp$age_effect_size)
    correlation_df = rbind(correlation_df,data.frame(matrix(c(sc,"DEG",test$estimate,test$p.value),nrow=1)))
}
colnames(correlation_df) = c("subclass","group","correlation","p-value")
correlation_df$correlation = as.numeric(correlation_df$correlation)
correlation_df$`p-value` = as.numeric(correlation_df$`p-value`)
correlation_df


write.csv(correlation_df,file=paste0(outdir,"Zeng_DEG_log2FC_and_Joint_Cabernet_genebody_5hmCG_diff_correlation.250423.csv"),quote=F,row.names=F)

var = "genebody"
mc = "CG"
dt = "5hmC"
if(dt=="true_5mC"){xlab=paste0("5mC_",var,"_",mc)}
if(dt=="5mC"){xlab=paste0("5mC+5hmC_",var,"_",mc)}
if(dt=="5hmC"){xlab=paste0("5hmC_",var,"_",mc)}
plots<-list()
g=1
for(sc in freq.df$Var1){
    print(paste0(dt,"_",mc,"_",var,"-",sc))
    pt<-plot_data[plot_data$subclass_label==sc,]
    pt_top = combined_plot_data[combined_plot_data$subclass_label==sc,]
    pt$group = "Down"
    pt$group[which(pt$age_effect_size > 0)] = "Up"
    
    x_range <- range(-abs(pt$diff),abs(pt$diff))
    y_range <- range(-abs(pt$age_effect_size),abs(pt$age_effect_size))
 
    x_margin <- (max(x_range) - min(x_range)) * 0.1  
    y_margin <- (max(y_range) - min(y_range)) * 0.1
 
    xlim <- c(min(x_range) - x_margin, max(x_range) + x_margin)
    ylim <- c(min(y_range) - y_margin, max(y_range) + y_margin)
    p1=ggplot(pt, aes(x=diff, y=age_effect_size)) +
                    geom_point(size = 1,color='#171717')+
                    theme_bw() +
                    scale_x_continuous(limits=xlim)+
                    scale_y_continuous(limits=ylim)+
                    theme(legend.position="none",
                    plot.title = element_text(hjust = 0.5,family="ArialMT"),
                    title = element_text(size=8,hjust = 0.5,family="ArialMT"),
                    axis.title.x = element_text(color="black", size=8,family="ArialMT"),
                    axis.title.y = element_text(color="black", size=8,family="ArialMT"),
                    axis.text.x = element_text(color="black", size=6,family="ArialMT"),
                    axis.text.y = element_text(color="black", size=6,family="ArialMT"),
                    axis.line=element_blank(),
                    axis.ticks=element_blank(),
                    panel.border=element_blank(),
                    panel.grid=element_blank(),
                    element_line(linetype = "dashed"), 
                    panel.grid.minor = element_blank())+  
                    labs(title=paste0(sc,"\nR: ",round(correlation_df$correlation[g],2),"\n","p-value: ",round(correlation_df$`p-value`[g],4)))+
                    xlab("genebody 5hmCG diff") +
                    ylab("RNA log2FC")+
                    geom_hline(yintercept=0, color = "black")+
                    geom_vline(xintercept=0, color = "black")+
                    geom_smooth(method = "lm", se = TRUE, color = "red") 
    plots[[g]]<-p1
    g=g+1
}

pdf("../../output/03.Aging_Mouse/07-aged_DEG/Zeng_DEG_RNA_log2FC_and_Joint_Cabernet_diff.color_by_age_effect_size.contain_NN.remove_dup_gene.without_axis_limits.20250421.pdf",height = 7,width = 12)
plot_grid(plotlist=plots[1:length(unique(freq.df$Var1))],ncol=6)
dev.off()


inte.col = readRDS("../../04.data/04.config_files/color.subclass.rds")

p3<-ggplot(correlation_df,aes(x=subclass,y=correlation,fill=subclass))+
    geom_bar(stat="identity",position = "dodge")+
    scale_x_discrete("subclass")+
    ylab("correlation")+
    scale_fill_manual(values = inte.col)+
    theme_bw()+
    scale_y_continuous(limits=c(0, 1), breaks=c(-1,-0.5,0,0.5,1),labels = function(x) sprintf("%.1f", x))+
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
pdf("../../output/03.Aging_Mouse/07-aged_DEG/Zeng_DEG_correlation_barplot.without_p_values.contain_NN.pdf",width=6,height=3)
print(p3)
dev.off()

correlation_df$subclass = factor(correlation_df$subclass,levels = subclass_order)
correlation_df$label <- ifelse(
  correlation_df$`p-value` <= 0.001, "***",
  ifelse(correlation_df$`p-value` <= 0.01 & correlation_df$`p-value` > 0.001, "**",
         ifelse(correlation_df$`p-value` <= 0.05 & correlation_df$`p-value` > 0.01, "*", "ns"))
)
p3<-ggplot(correlation_df,aes(x=subclass,y=correlation,fill=subclass))+
    geom_bar(stat="identity",position = "dodge")+
    scale_x_discrete("subclass")+
    ylab("correlation")+
    geom_text(data = correlation_df, aes(x = subclass, y = 0.8, label = correlation_df$label), 
            hjust = 0.5, vjust = -0.6, size = 2, color = "black") + # round(correlation_df$`p-value`,digit=2)
    scale_fill_manual(values = inte.col)+
    theme_bw()+
    scale_y_continuous(limits=c(0, 1), breaks=c(-1,-0.5,0,0.5,1),labels = function(x) sprintf("%.1f", x))+
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

pdf("../../output/03.Aging_Mouse/07-aged_DEG/Zeng_DEG_correlation_barplot.with_p_values.contain_NN.remove_dup_gene.20250417.pdf",width=6,height=3)
print(p3)
dev.off()