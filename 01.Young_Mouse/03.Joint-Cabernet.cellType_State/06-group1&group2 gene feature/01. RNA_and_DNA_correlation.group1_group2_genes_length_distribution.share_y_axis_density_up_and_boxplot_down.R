##### 01. import packages #####
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

##### 02. data process #####
group1 = read.table("group1_gene_list.txt",header=F)
head(group1)

mm10_df = fread("../../../03.data/01.ref/mm10.genes.bed",header=T,sep="\t",data.table=F)
mm10_df[1:3,]

colnames(group1) = "genename"
group1$geneid = mm10_df[match(group1$genename,mm10_df$gene_name),"gene_id"]
head(group1)

group1$start = mm10_df[match(group1$genename,mm10_df$gene_name),"start"]
group1$end = mm10_df[match(group1$genename,mm10_df$gene_name),"end"]
group1$gene_length = group1$end - group1$start


CpG_number = read.csv("../../../03.data/01.ref/gene_CpG_number_metainfo.csv",header=T)
CpG_number[1:5,]

group1$Cpg_number = CpG_number[match(group1$geneid,CpG_number$gene_id),"Cpg_number"]
head(group1)

group1$log_gene_length = log10(group1$gene_length)

group2 = read.table("group2_gene_list.txt",header=F)
colnames(group2) = "genename"
group2$geneid = mm10_df[match(group2$genename,mm10_df$gene_name),"gene_id"]
group2$start = mm10_df[match(group2$genename,mm10_df$gene_name),"start"]
group2$end = mm10_df[match(group2$genename,mm10_df$gene_name),"end"]
group2$gene_length = group2$end - group2$start
group2$Cpg_number = CpG_number[match(group2$geneid,CpG_number$gene_id),"Cpg_number"]
group2$log_gene_length = log10(group2$gene_length)
head(group2)

group1$group = "Group1"
group2$group = "Group2"
combine_df = rbind(group1,group2)

combine_df$log_CpG = log10(combine_df$Cpg_number)

inte.col = c("Group1"="#faada7","Group2"="#65d9dc")
color.col = c("Group1"="#ef7773","Group2"="#02bfbf")

##### 03. plot #####
##### gene length #####
pdf("log_gene_length.up_and_down.pdf",width=7,height=4)
ks_test=ks.test(combine_df[combine_df$group == "Group1","log_gene_length"],combine_df[combine_df$group == "Group2","log_gene_length"])
pval =ks_test$p.value
density_plot <- ggplot(combine_df, aes(x=log_gene_length))+
    geom_density(aes(fill=group,color=group), alpha=0.6)+ 
    xlab('log10 gene length')+
    ylab('Density')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(), 
            axis.text.y = element_text(size=14,face="bold",color = "black"),
            legend.key.size = unit(20, "pt"),
            legend.title = element_text(face="bold",size=15),
            legend.text = element_text(face="bold",size=15),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=18))+
    scale_x_continuous(limits=c(4, 6.5),breaks=c(4,4.5,5,5.5,6)) 
boxplot_plot <- ggplot(combine_df, aes(x = group, y = log_gene_length, fill = group)) +
  geom_boxplot(width = 0.4) +
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(), 
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.text.y = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14),
        legend.position = "none")+
  scale_y_continuous(limits=c(4, 6.5),breaks=c(4,4.5,5,5.5,6))+coord_flip()

# patchwork to combine plots
combined_plot <- density_plot / boxplot_plot + 
  plot_layout(ncol = 1,heights = c(2, 1))
# add title
final_plot <- combined_plot +
  plot_annotation(title ="Gene length",theme = theme(plot.title = element_text(hjust=0.5,face="bold", size=18)))

print(final_plot)
dev.off()

head(combine_df)

##### genebody_CpG_number #####
pdf("log_genebody_CpG_number.up_and_down.pdf",width=7,height=4)
ks_test=ks.test(combine_df[combine_df$group == "Group1","log_CpG"],combine_df[combine_df$group == "Group2","log_CpG"])
pval =ks_test$p.value
density_plot <- ggplot(combine_df, aes(x=log_CpG))+
    geom_density(aes(fill=group,color=group), alpha=0.6)+ 
    xlab('log10 genebody CpG number')+
    ylab('Density')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),                 
            axis.text.x = element_blank(),  
            axis.text.y = element_text(size=14,face="bold",color = "black"),
            legend.key.size = unit(20, "pt"),
            legend.title = element_text(face="bold",size=15),
            legend.text = element_text(face="bold",size=15),
            text = element_text(face="bold",size = 10),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face="bold", size=18))+
    scale_x_continuous(limits=c(2, 4.5),breaks=c(2,2.5,3,3.5,4,4.5)) 
boxplot_plot <- ggplot(combine_df, aes(x = group, y = log_CpG, fill = group)) +
  geom_boxplot(width = 0.4) +
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(), 
        axis.line = element_line(linewidth=0.5),
        axis.title = element_blank(),
        axis.text = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14),
        legend.position = "none")+
  scale_y_continuous(limits=c(2, 4.5),breaks=c(2,2.5,3,3.5,4,4.5))+coord_flip()

# patchwork to combine plots
combined_plot <- density_plot / boxplot_plot + 
  plot_layout(ncol = 1,heights = c(2, 1))
# add title
final_plot <- combined_plot +
    plot_annotation(title ="CpG number",theme = theme(plot.title = element_text(hjust=0.5,face="bold", size=18)))

print(final_plot)
dev.off()


