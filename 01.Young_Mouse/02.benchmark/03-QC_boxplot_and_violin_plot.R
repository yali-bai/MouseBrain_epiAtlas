library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)

DNA_stat = read.csv("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.DNA_QC_stat.young.csv",header=T)
DNA_stat = DNA_stat[DNA_stat$total_QC == 1,]
DNA_stat$group = "5mC+5hmC"
DNA_stat$group[which(DNA_stat$Library == "hmC")] = "5hmC"
DNA_stat$group = factor(DNA_stat$group,levels = c("5hmC","5mC+5hmC"))

DNA_stat$COVERAGE = as.numeric(DNA_stat$COVERAGE)*100

DNA_stat_no_outliers <- DNA_stat %>%
  group_by(group) %>%
  mutate(
    Q1 = quantile(COVERAGE, 0.25, na.rm = TRUE),
    Q3 = quantile(COVERAGE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(COVERAGE >= lower_bound & COVERAGE <= upper_bound) %>%
  ungroup()

# 
summary_df_no_outliers <- DNA_stat_no_outliers %>%
  group_by(group) %>%
  summarise(
    mean = mean(COVERAGE, na.rm = TRUE),
    sd = sd(COVERAGE, na.rm = TRUE),
    se = sd(COVERAGE, na.rm = TRUE) / sqrt(n()),
    median = median(COVERAGE, na.rm = TRUE),
    Q1 = quantile(COVERAGE, 0.25, na.rm = TRUE),
    Q3 = quantile(COVERAGE, 0.75, na.rm = TRUE),
    n = n(),
    max = max(COVERAGE),
    min = min(COVERAGE)
  )


pdf("genome_coverage.boxplot.pdf",width = 4,height = 4)
ggplot() +
    geom_boxplot(data = DNA_stat,aes(x = group, y = COVERAGE,fill = group),linewidth=0.5,na.rm = TRUE,outliers = FALSE)+
    geom_errorbar(data = summary_df_no_outliers,aes(ymin = min, ymax = max, x = group),width = 0.2,color = "black",linewidth = 0.5) +
    scale_fill_manual(values = c('5hmC'='#e47927',"5mC+5hmC"='#1471aa'))+
    labs(y = "Genome coverage%") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text = element_text(size=20,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=25,face="bold",color = "black"))+
    scale_y_continuous(limits=c(0,20),breaks=seq(0,20,5))
dev.off()


DNA_stat$dna_CGn = as.numeric(DNA_stat$dna_CGn)/1e6
DNA_stat_no_outliers <- DNA_stat %>%
  group_by(group) %>%
  mutate(
    Q1 = quantile(dna_CGn, 0.25, na.rm = TRUE),
    Q3 = quantile(dna_CGn, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(dna_CGn >= lower_bound & dna_CGn <= upper_bound) %>%
  ungroup()

#
summary_df_no_outliers <- DNA_stat_no_outliers %>%
  group_by(group) %>%
  summarise(
    mean = mean(dna_CGn, na.rm = TRUE),
    sd = sd(dna_CGn, na.rm = TRUE),
    se = sd(dna_CGn, na.rm = TRUE) / sqrt(n()),
    median = median(dna_CGn, na.rm = TRUE),
    Q1 = quantile(dna_CGn, 0.25, na.rm = TRUE),
    Q3 = quantile(dna_CGn, 0.75, na.rm = TRUE),
    n = n(),
    max = max(dna_CGn),
    min = min(dna_CGn)
  )
pdf("detect_CG_sites.boxplot.pdf",width = 4,height = 4.5)
ggplot() +
    geom_boxplot(data = DNA_stat,aes(x = group, y = dna_CGn,fill = group),linewidth=0.5,na.rm = TRUE,outliers = FALSE)+
    geom_errorbar(data = summary_df_no_outliers,aes(ymin = min, ymax = max, x = group),width = 0.2,color = "black",linewidth = 0.5) +
    scale_fill_manual(values = c('5hmC'='#e47927',"5mC+5hmC"='#1471aa'))+
    labs(y = "Detected CG sites") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text = element_text(size=20,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=25,face="bold",color = "black"),
            plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "mm"))+
    scale_y_continuous(limits=c(0,4),breaks=c(0,1,2,3,4),labels = c(0,1,2,3,4))+
    coord_cartesian(clip = "off") +
    annotate("text", x = -Inf, y = Inf, 
           label = "x10^6", 
           hjust = 0.5, vjust = -0.5,  
           size = 7, fontface = "bold", parse = TRUE)
dev.off()



RNA_stat = read.csv("../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.RNA_QC_stat.young.csv",header=T)
RNA_stat = RNA_stat[RNA_stat$QC_after_integration ==1,]

RNA_stat$Gene_gene_number = as.numeric(RNA_stat$Gene_gene_number)/1e3

pdf("Gene_number.boxplot.pdf",width = 4,height = 4.5)
ggplot(data = RNA_stat,aes(x = Neuron_non_neuron, y = Gene_gene_number,fill = Neuron_non_neuron)) +
    geom_violin(linewidth=0.2,width=0.8,na.rm = TRUE)+
    geom_boxplot(linewidth=0.2,width=0.15,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = c('Neuron'='#628ebf',"Non_neuron"='#f1bc24'))+
    labs(y = "Number of genes") +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.text = element_text(size=20,face="bold",color = "black"),
            legend.position = "none",
            text = element_text(face="bold",size = 10),
            plot.title = element_text(size = 18,face="bold",hjust = 0.5),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.line.y = element_line(linewidth=0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=25,face="bold",color = "black"),
            plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "mm"))+
    coord_cartesian(clip = "off") +
    annotate("text", x = -Inf, y = Inf, 
           label = "x10^3", 
           hjust = 0.5, vjust = -0.5,  
           size = 7, fontface = "bold", parse = TRUE)
dev.off()

gene_number = read.csv("RNA_gene_number_in_different_sequencing_technology.csv",header=T)
Joint_Cabernet_RNA_stat = read.csv("TSO-joint.RNA_QC_stat.young.add_celltype.csv",header=T)
Joint_Cabernet_RNA_stat = Joint_Cabernet_RNA_stat[Joint_Cabernet_RNA_stat$QC_after_integration == 1,]

plot_df = data.frame(gene_number = Joint_Cabernet_RNA_stat$Gene_gene_number, tech = "Joint_Cabernet")
plot_df = rbind(plot_df,data.frame(gene_number = gene_number[,3], tech = "10x"))
plot_df = rbind(plot_df,data.frame(gene_number = gene_number[,4], tech = "Smart"))
plot_df$tech = factor(plot_df$tech,levels = c("Smart","10x","Joint_Cabernet"))

DNA_stat_no_outliers <- plot_df %>%
  group_by(tech) %>%
  mutate(
    Q1 = quantile(gene_number, 0.25, na.rm = TRUE),
    Q3 = quantile(gene_number, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(gene_number >= lower_bound & gene_number <= upper_bound) %>%
  ungroup()

# 
summary_df_no_outliers <- DNA_stat_no_outliers %>%
  group_by(tech) %>%
  summarise(
    mean = mean(gene_number, na.rm = TRUE),
    sd = sd(gene_number, na.rm = TRUE),
    se = sd(gene_number, na.rm = TRUE) / sqrt(n()),
    median = median(gene_number, na.rm = TRUE),
    Q1 = quantile(gene_number, 0.25, na.rm = TRUE),
    Q3 = quantile(gene_number, 0.75, na.rm = TRUE),
    n = n(),
    max = max(gene_number),
    min = min(gene_number)
  )


pdf("Gene_number_of_different_sequencing_tech.boxplot.pdf",width = 6.5,height = 3)
ggplot() +
    geom_errorbar(data = summary_df_no_outliers,aes(ymin = min, ymax = max, x = tech),width = 0.2,color = "black",linewidth = 0.5)+
    geom_boxplot(data = plot_df,aes(x = tech, y = gene_number,fill = "#cccccc"),linewidth=0.5,width=0.5,na.rm = TRUE,outliers = FALSE)+
    scale_fill_manual(values = "#cccccc")+
    scale_color_manual(values = "black")+
    labs(y = "Number of genes") +
    scale_y_continuous(breaks = seq(2000,14000,2000))+
    scale_x_discrete(labels = c(paste0("Smart-seq\n(Yao et al .)"),paste0("10x RNA-seq\n(Yao et al .)"),"Joint_Cabernet"))+
    theme_bw()+
    theme(panel.grid.major.y=element_line(color = "grey90", linewidth = 0.3), 
            panel.grid.major.x=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            axis.text.y = element_text(size=15,color = "black"),
            axis.text.x = element_text(angle=60,vjust = 1,hjust =1,size=15,color = "black"),
            legend.position = "none",
            text = element_text(size = 15),
            plot.title = element_text(size = 20,hjust = 0.5),
            axis.ticks.x = element_line(linewidth=0.5),
            axis.line = element_line(linewidth=0.5),
            axis.title.x = element_text(size=20,color = "black"),
            axis.title.y = element_blank())+
    coord_flip()#+
    
dev.off()




