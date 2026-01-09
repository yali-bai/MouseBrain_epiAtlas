##### 01. import packages #####
library(ggpubr)
library(Seurat)
library(dplyr)
library(Matrix)
library(dplyr)
library(reshape2)
library(tidyverse)

##### 02. data process #####
all_gene_set <- read.csv("all_gene_upset_group_list.csv")
head(all_gene_set);dim(all_gene_set)

unique(all_gene_set$group)
unique(all_gene_set$gene_type)

group1_genes <- all_gene_set[all_gene_set$group == '[5mCG -]:[5hmCG +]:[(5mCG+5hmCG) -]',]$gene_id
group1_genes;length(group1_genes);length(unique(group1_genes))

group2_genes <- all_gene_set[all_gene_set$group == '[5mCG -]:[5hmCG +]',]$gene_id
group2_genes;length(group2_genes);length(unique(group2_genes))

##### subclass RNA calculation #####
seuratObj <- readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/Joint_Cabernet.with_celltype.rds")
head(seuratObj@meta.data)

DefaultAssay(seuratObj) = "RNA"
seuratObj = NormalizeData(seuratObj)
our_dat <- seuratObj@assays$RNA$data

our_label <- seuratObj@meta.data %>% dplyr::select(major_class.annotated_all,subclass.annotated_all) %>% mutate(cell=rownames(.))
head(our_label)

subclass <- readRDS("../../../03.data/04.config_files/subclass_order.rds")

subclass = intersect(subclass,unique(our_label$subclass.annotated_all))

our_label <- our_label[our_label$subclass.annotated_all %in% subclass,]

our_dat <- our_dat[,rownames(our_label)]
our_dat_t <- t(our_dat)
our_subclass <- our_label[match(rownames(our_dat_t),rownames(our_label)),"subclass.annotated_all"]
our_dat_final <- aggregate(our_dat_t, by=list(our_subclass), FUN = function(x) mean(x, na.rm = TRUE))
our_dat_final <- column_to_rownames(our_dat_final,'Group.1')
our_dat_final <- t(our_dat_final)
head(our_dat_final)
dim(our_dat_final)


saveRDS(our_dat_final,"02-our_subclass_mean_dat_final.rds")
write.csv(our_dat_final, file = "02-our_subclass_mean_dat_final.csv", quote=F)


young_RNA <- our_dat_final
head(young_RNA);dim(young_RNA)

table(group1_genes %in% rownames(young_RNA));table(group2_genes %in% rownames(young_RNA))

library(ggplot2)
library(tidyr)
library(patchwork)

group1_expr <- young_RNA[group1_genes, ]
group2_expr <- young_RNA[group2_genes, ]
group1_sd <- apply(group1_expr, 1, sd)
group2_sd <- apply(group2_expr, 1, sd)
group1_cv <- apply(group1_expr, 1, function(x) sd(x) / mean(x))
group2_cv <- apply(group2_expr, 1, function(x) sd(x) / mean(x))
    

sd_df <- data.frame(
  value = c(group1_sd, group2_sd),
  metric = "SD",
  group = rep(c("group1", "group2"), c(length(group1_sd), length(group2_sd)))
)

cv_df <- data.frame(
  value = c(group1_cv, group2_cv),
  metric = "CV",
  group = rep(c("group1", "group2"), c(length(group1_cv), length(group2_cv)))
)

plot_df <- bind_rows(sd_df, cv_df)
plot_df;dim(plot_df)

##### 03. plot #####
pdf("group1_group2_genes_RNA_expression_across_subclasses_CV_SD.pdf",width = 6,height = 4)
ggplot(plot_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, color = "black") +
  facet_wrap(~ metric, scales = "free_y") +
  #scale_fill_brewer(palette = "Set2") +
  ylab("Expression Variability") +
  xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    comparisons = list(c("group1", "group2")),
    tip.length = 0.01,
    size = 4
  )
dev.off()
