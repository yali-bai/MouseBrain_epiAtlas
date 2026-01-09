library(ggplot2)
library(cowplot)
library(Seurat)
library(dplyr)
library(future)
library(presto)
library(stringr)
library(getopt)
library(data.table)
library(ggunchained)
library(reshape2)
library(ggpubr)


setwd("./all_cell_correlation")

hmC = read.csv("all_cells_5hmC_CG_genebody_gene_correlation_results.csv",header=T)
hmC[1:3,]

mC = read.csv("all_cells_5mC_CG_genebody_gene_correlation_results.csv",header=T)
true_mC = read.csv("all_cells_true_5mC_CG_genebody_gene_correlation_results.csv",header=T)

merge_df = rbind(hmC,mC)
merge_df = rbind(merge_df,true_mC)

merge_df$group = NA
merge_df$group[which(merge_df$datatype == "5hmC")] = "5hmCG"
merge_df$group[which(merge_df$datatype == "5mC")] = "5mCG+5hmCG"
merge_df$group[which(merge_df$datatype == "true_5mC")] = "5mCG"
unique(merge_df$group)
head(merge_df)

unique(merge_df$datatype)

merge_df$group = factor(merge_df$group,levels=c('5mCG+5hmCG','5mCG','5hmCG'))

p<-ggplot(merge_df, aes(x = Correlation, group = group, color = group,fill = group,alpha=0.9)) +  
                geom_density()+
                scale_fill_manual(values = c("5hmCG" = "#d07431", "5mCG" = "#3076c3", "5mCG+5hmCG" = "#777a81"))+# Set the transparency of each group separately 
                scale_color_manual(values = c("5hmCG" = "#d07431", "5mCG" = "#3076c3", "5mCG+5hmCG" = "#777a81"))+
                scale_x_continuous(breaks=c(-0.2,0,0.2),limit=c(-0.2,0.2))+
                scale_y_continuous(breaks = seq(from = 0, to = 20, by = 10),  
                     limits = c(0, 20))+             
                labs(x ="PCC between gene body methylation and expression" ,y="Density of genes")+
                theme_bw()+
                theme(legend.position= "right",
                plot.title = element_text(hjust = 0.5,family="ArialMT"),
                title = element_text(size=7,hjust = 0.5,family="ArialMT"),
                axis.title.x = element_text(color="black", size=19,family="ArialMT"),
                axis.title.y = element_text(color="black", size=25,family="ArialMT"),
                axis.text.x = element_text(color="black", size=20,family="ArialMT", hjust = 1,),
                axis.text.y = element_text(color="black", size=20,family="ArialMT"),
                legend.title = element_text(color="black", size=9,family="ArialMT"),
                axis.line=element_line(linewidth=0.6),
                axis.ticks=element_line(linewidth=0.6),
                panel.border=element_blank(),
                panel.grid=element_blank())

pdf("PCC_between_gene_body_methylation_and_expression.density_plot.pdf",width=7,height=4)
print(p)
dev.off()


p<-ggplot(merge_df, aes(x = Correlation, group = group, color = group,fill = group)) +  
                geom_density(linewidth = 1.5,alpha = 0.6)+
                scale_fill_manual(values = c("5hmCG" = "#d07431", "5mCG" = "#3076c3", "5mCG+5hmCG" = "#777a81"))+# Set the transparency of each group separately 
                scale_color_manual(values = c("5hmCG" = "#d07431", "5mCG" = "#3076c3", "5mCG+5hmCG" = "#777a81"))+
                scale_x_continuous(breaks=c(-0.2,0,0.2),limit=c(-0.2,0.2))+
                scale_y_continuous(breaks = seq(from = 0, to = 20, by = 10),  
                     limits = c(0, 20))+             
                labs(x ="PCC between gene body methylation and expression" ,y="Density of genes")+
                theme_bw()+
                theme(legend.position= "right",
                plot.title = element_text(hjust = 0.5,family="ArialMT"),
                title = element_text(size=7,hjust = 0.5,family="ArialMT"),
                axis.title.x = element_text(color="black", size=19,family="ArialMT"),
                axis.title.y = element_text(color="black", size=25,family="ArialMT"),
                axis.text.x = element_text(color="black", size=20,family="ArialMT", hjust = 1,),
                axis.text.y = element_text(color="black", size=20,family="ArialMT"),
                legend.title = element_blank(),
                legend.text = element_text(color="black", size=12,family="ArialMT"),
                axis.line=element_line(linewidth=0.6),
                axis.ticks=element_line(linewidth=0.6),
                panel.border=element_blank(),
                panel.grid=element_blank())+
                guides(fill = guide_legend(override.aes = list(alpha = 0.6)))

pdf("PCC_between_gene_body_methylation_and_expression.density_plot.pdf",width=8,height=4)
print(p)
dev.off()



