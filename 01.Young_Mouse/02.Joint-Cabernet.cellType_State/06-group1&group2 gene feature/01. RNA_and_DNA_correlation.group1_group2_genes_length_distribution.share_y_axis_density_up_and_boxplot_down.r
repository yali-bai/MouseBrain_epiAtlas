##### 01. import packages #####
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

##### 02. set working path #####
combine_df = read.csv("../../../04.data/05.intermediate_files/01.RNA/01.Young_Mouse/gene_metainfo.of_group1_group2.csv")
head(combine_df)

inte.col = c("Group1"="#f1a0a0","Group2"="#7cafd1")
color.col = c("Group1"="#e46363","Group2"="#3078ab")
pdf("../../../output/01.Young_Mouse/02-correlation_calculation/group1&group2 gene feature/log_gene_length.up_and_down.pdf",width=7,height=4)

ks_test=ks.test(combine_df[combine_df$group == "Group1","log_gene_length"],combine_df[combine_df$group == "Group2","log_gene_length"])
pval =ks_test$p.value
density_plot <- ggplot(combine_df, aes(x=log_gene_length))+
    geom_density(aes(fill=group,color=group), alpha=0.6)+ 
    xlab('log10 gene length')+
    scale_fill_manual(values = inte.col)+
    scale_color_manual(values = color.col)+
    theme_bw()+
    geom_vline(xintercept=mean(combine_df[combine_df$group == "Group1","log_gene_length"]), color = "#e46363", linetype = "dashed")+
    geom_vline(xintercept=mean(combine_df[combine_df$group == "Group2","log_gene_length"]), color = "#3078ab", linetype = "dashed")+
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
    scale_x_continuous(limits=c(4, 6.5),breaks=c(4,4.5,5,5.5,6)) 
boxplot_plot <- ggplot(combine_df, aes(x = group, y = log_gene_length, fill = group)) +
  geom_boxplot(width = 0.4) +
  theme_bw()+
  geom_hline(yintercept=mean(combine_df[combine_df$group == "Group1","log_gene_length"]), color = "#e46363", linetype = "dashed")+
  geom_hline(yintercept=mean(combine_df[combine_df$group == "Group2","log_gene_length"]), color = "#3078ab", linetype = "dashed")+
  scale_fill_manual(values = inte.col)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(), 
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.text.y = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14), #element_blank(),  # 去掉y轴文本
        axis.title = element_blank(),
        axis.text.x = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14),
        legend.position = "none")+
  scale_y_continuous(limits=c(4, 6.5),breaks=c(4,4.5,5,5.5,6))+coord_flip()

# patchwork
combined_plot <- density_plot / boxplot_plot + 
  plot_layout(ncol = 1,heights = c(2, 1))

final_plot <- combined_plot +
  plot_annotation(title =paste0("KS-test p-value: ",round(pval,digit=3)),theme = theme(plot.title = element_text(hjust=0.5,face="bold", size=18)))

print(final_plot)
dev.off()

pdf("../../../output/01.Young_Mouse/02-correlation_calculation/group1&group2 gene feature/log_genebody_CpG_number.up_and_down.pdf",width=7,height=4)

ks_test=ks.test(combine_df[combine_df$group == "Group1","log_CpG"],combine_df[combine_df$group == "Group2","log_CpG"])
pval =ks_test$p.value
density_plot <- ggplot(combine_df, aes(x=log_CpG))+
    geom_density(aes(fill=group,color=group), alpha=0.6)+ 
    xlab('log10 genebody CpG number')+
    scale_fill_manual(values = inte.col)+
    scale_color_manual(values = color.col)+
    theme_bw()+
    geom_vline(xintercept=mean(combine_df[combine_df$group == "Group1","log_CpG"]), color = "#e46363", linetype = "dashed")+
    geom_vline(xintercept=mean(combine_df[combine_df$group == "Group2","log_CpG"]), color = "#3078ab", linetype = "dashed")+
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
  geom_hline(yintercept=mean(combine_df[combine_df$group == "Group1","log_CpG"]), color = "#e46363", linetype = "dashed")+
  geom_hline(yintercept=mean(combine_df[combine_df$group == "Group2","log_CpG"]), color = "#3078ab", linetype = "dashed")+
  scale_fill_manual(values = inte.col)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(), 
        axis.line = element_line(linewidth=0.5),
        axis.title = element_blank(),
        axis.text = element_text(angle=0,vjust = 1,hjust =0.5,color = "black",face="bold",size=14),
        legend.position = "none")+
  scale_y_continuous(limits=c(2, 4.5),breaks=c(2,2.5,3,3.5,4,4.5))+coord_flip()

# patchwork
combined_plot <- density_plot / boxplot_plot + 
  plot_layout(ncol = 1,heights = c(2, 1))

final_plot <- combined_plot +
  plot_annotation(title =paste0("KS-test p-value: ",round(pval,digit=3)),theme = theme(plot.title = element_text(hjust=0.5,face="bold", size=18)))


print(final_plot)
dev.off()


