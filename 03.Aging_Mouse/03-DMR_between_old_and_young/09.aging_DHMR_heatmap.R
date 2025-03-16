#########    All "our" in the following code refers to Joint Cabernet.
##### 01.import packages #####
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library","/share/analysisdata/Methyl/public/rna/lib/R/library"))
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(ggpointdensity) 
library(cowplot)
library(ggtext)
library(ggpubr)
library(ggunchained)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
# setwd("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/05.significant_DMR_DHMR_mcds/20250217_update.filter_before_calculating_p_value_adjust")

##### 03.set legend color #####
#inte.col = readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/color.subclass.rds")
inte.col = readRDS("../../input/03-aging/color.subclass.rds")

##### 04.read subclass order #####
#subclass_order = readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/order.subclass.rds")
subclass_order = readRDS("../../input/03-aging/order.subclass.rds")

##### 05.plot #####
## read input file which contains 5hmCG, 5mCG, 5mCG_5hmCG mean diff value of all segments in all subclasses ##
merge.df.dcast = readRDS("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/05.significant_DMR_DHMR_mcds/20250217_update.filter_before_calculating_p_value_adjust/merge.df.dcast.rds")
head(merge.df.dcast)

## significant DHMRs ##
#DHMR_result=readRDS("../../04.DMR_DHMR_calculate_and_filter/20250217_update.filter_before_calculating_p_value_adjust/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")
DHMR_result=readRDS("../../output/03-aging/03-DMRs_DHMRs/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")
head(DHMR_result)

## select hyper DHMRs and sort by subclass order and diff value in descending order ##
DHMR_result.sig = subset(DHMR_result, diff > 0.05 )
DHMR_result.sig$cluster = factor(DHMR_result.sig$cluster,levels = intersect(subclass_order,unique(DHMR_result.sig$cluster)))
DHMR_result.sig.sorted <- DHMR_result.sig[order(DHMR_result.sig$cluster, DHMR_result.sig$diff), ]
DHMR_result.sig.sorted[1:5,]

## remove DHMRs which appear in more than one subclass, that is , retain subclass-specific DHMRs ##
freq.df = as.data.frame(table(DHMR_result.sig.sorted$chrom))
head(freq.df)
repeat.v = freq.df[which(freq.df$Freq > 1),"Var1"]
remove.idx = which(DHMR_result.sig.sorted$chrom %in% repeat.v) # index of DHMRs which appear in more than one subclass
DHMR_result.sig.sorted.rm_repeat = DHMR_result.sig.sorted[-remove.idx,]

## extract top1000 subclass-specific DHMRs ##
top_df <- DHMR_result.sig.sorted.rm_repeat %>%
  group_by(cluster) %>%
  arrange(desc(diff)) %>%
  slice(1:1000)
head(top_df)

## generate data for plotting heatmap of 5hmCG change ##
total_mC.dcast = merge.df.dcast[which(merge.df.dcast$type == "5hmCG"),]
total_mC.dcast_subset = total_mC.dcast[,c("segment","subclass","diff")]
total_mC.diff_matrix = dcast(total_mC.dcast_subset,segment~subclass)
rownames(total_mC.diff_matrix) = total_mC.diff_matrix$segment
total_mC.diff_matrix = total_mC.diff_matrix[,-1]
total_mC.diff_matrix[1:5,1:5]
total_mC.diff_matrix.sorted = total_mC.diff_matrix[top_df$chrom,intersect(subclass_order,unique(merge.df.dcast$subclass))]
total_mC.diff_matrix.sorted[1:5,1:5]
row_order = rownames(total_mC.diff_matrix.sorted)
col_order = colnames(total_mC.diff_matrix.sorted)
write.csv(total_mC.diff_matrix.sorted,file = "../../output/03-aging/DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5hmCG_mean_diff_old_minus_young_matrix.subclass_specific_top_500.csv",quote=F,row.names=T,col.names=T)

## generate data for plotting heatmap of 5mCG_5hmCG change ##
total_mC.dcast = merge.df.dcast[which(merge.df.dcast$type == "5mCG_5hmCG"),]
total_mC.dcast_subset = total_mC.dcast[,c("segment","subclass","diff")]
total_mC.diff_matrix = dcast(total_mC.dcast_subset,segment~subclass)
rownames(total_mC.diff_matrix) = total_mC.diff_matrix$segment
total_mC.diff_matrix = total_mC.diff_matrix[,-1]
total_mC.diff_matrix.sorted = total_mC.diff_matrix[row_order,col_order]
write.csv(total_mC.diff_matrix.sorted,file = "../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5mCG_5hmCG_mean_diff_old_minus_young_matrix.subclass_specific_top_500.csv",quote=F,row.names=T,col.names=T)

## generate data for plotting heatmap of 5mCG change ##
total_mC.dcast = merge.df.dcast[which(merge.df.dcast$type == "5mCG"),]
total_mC.dcast_subset = total_mC.dcast[,c("segment","subclass","diff")]
total_mC.diff_matrix = dcast(total_mC.dcast_subset,segment~subclass)
rownames(total_mC.diff_matrix) = total_mC.diff_matrix$segment
total_mC.diff_matrix = total_mC.diff_matrix[,-1]
total_mC.diff_matrix.sorted = total_mC.diff_matrix[row_order,col_order]
write.csv(total_mC.diff_matrix.sorted,file = "../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5mCG_mean_diff_old_minus_young_matrix.subclass_specific_top_500.csv",quote=F,row.names=T,col.names=T)

## plot 5mCG mean diff value heatmap ##
merge.df.dcast.temp = fread("../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5mCG_mean_diff_old_minus_young_matrix.subclass_specific_top_1000.csv",data.table=F,header=T)
rownames(merge.df.dcast.temp) = merge.df.dcast.temp[,1]
merge.df.dcast.temp  = merge.df.dcast.temp[,-1]
total_mC.diff_matrix = as.matrix(merge.df.dcast.temp)
ha1=HeatmapAnnotation(subclass = colnames(total_mC.diff_matrix),show_legend = FALSE,annotation_name_side = "right",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))

#pdf("../../result/aging/DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.true_5mCG_diff_old_minus_young.heatmap.pdf",width=8.2677,height=21.2598)
png(filename = "../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.true_5mCG_diff_old_minus_young.heatmap.png", width = 8.2677, height = 21.2598, units = "cm", res = 300)
col_fun =colorRamp2(c(-0.2,-0.1,0,0.1,0.2), c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3")) 
Heatmap(total_mC.diff_matrix,
        column_gap = unit(0, "points"),
        bottom_annotation = ha1,
        col = col_fun,
        cluster_rows=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names=F,
        show_row_names=F,
        row_names_side="left",
        column_names_side="bottom",
        column_names_rot = 75,
        column_dend_side = "top",
        border = TRUE,
        row_title = "Aging hyper-DHMRs (Aged > Young)",
        column_title = "5mCG change",
        column_title_gp = gpar(fontsize = 25, fontface = "bold"),
        column_title_side = "top",
        column_title_rot = 0,
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 25, fontface = "bold"),
        show_heatmap_legend=FALSE,
        use_raster = FALSE)
dev.off()

## plot 5mCG_5hmCG mean diff value heatmap ##
merge.df.dcast.temp = fread("../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5mCG_5hmCG_mean_diff_old_minus_young_matrix.subclass_specific_top_1000.csv",data.table=F,header=T)
rownames(merge.df.dcast.temp) = merge.df.dcast.temp[,1]
merge.df.dcast.temp  = merge.df.dcast.temp[,-1]
total_mC.diff_matrix = as.matrix(merge.df.dcast.temp)
ha1=HeatmapAnnotation(subclass = colnames(total_mC.diff_matrix),show_legend = FALSE,annotation_name_side = "right",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))

#pdf("../../result/aging/DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.5mCG_5hmCG_diff_old_minus_young.heatmap.pdf",width=8.2677,height=21.2598)
png(filename = "../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.5mCG_5hmCG_diff_old_minus_young.heatmap.png", width = 8.2677, height = 21.2598, units = "cm", res = 300)
col_fun =colorRamp2(c(-0.2,-0.1,0,0.1,0.2), c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3")) 
Heatmap(total_mC.diff_matrix,
        column_gap = unit(0, "points"),
        bottom_annotation = ha1,
        col = col_fun,
        cluster_rows=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names=F,
        show_row_names=F,
        row_names_side="left",
        column_names_side="bottom",
        column_names_rot = 75,
        column_dend_side = "top",
        border = TRUE,
        row_title = "Aging hyper-DHMRs (Aged > Young)",
        column_title = "5mCG_5hmCG change",
        column_title_gp = gpar(fontsize = 25, fontface = "bold"),
        column_title_side = "top",
        column_title_rot = 0,
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 25, fontface = "bold"),
        show_heatmap_legend=FALSE,
        use_raster = FALSE)
dev.off()

## plot 5hmCG mean diff value heatmap ##
merge.df.dcast.temp = fread("../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/hyper_DHMR.5hmCG_mean_diff_old_minus_young_matrix.subclass_specific_top_1000.csv",data.table=F,header=T)
rownames(merge.df.dcast.temp) = merge.df.dcast.temp[,1]
merge.df.dcast.temp  = merge.df.dcast.temp[,-1]
total_mC.diff_matrix = as.matrix(merge.df.dcast.temp)
ha1=HeatmapAnnotation(subclass = colnames(total_mC.diff_matrix),show_legend = FALSE,annotation_name_side = "right",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))

#pdf("../../result/aging/DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.5hmCG_diff_old_minus_young.heatmap.pdf",width=8.2677,height=21.2598)
png(filename = "../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/aging_hyper_DHMRs_top1000.5hmCG_diff_old_minus_young.heatmap.png", width = 8.2677, height = 21.2598, units = "cm", res = 300)
col_fun =colorRamp2(c(-0.2,-0.1,0,0.1,0.2), c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3")) 
Heatmap(total_mC.diff_matrix,
        column_gap = unit(0, "points"),
        bottom_annotation = ha1,
        col = col_fun,
        cluster_rows=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names=F,
        show_row_names=F,
        row_names_side="left",
        column_names_side="bottom",
        column_names_rot = 75,
        column_dend_side = "top",
        border = TRUE,
        row_title = "Aging hyper-DHMRs (Aged > Young)",
        column_title = "5hmCG change",
        column_title_gp = gpar(fontsize = 25, fontface = "bold"),
        column_title_side = "top",
        column_title_rot = 0,
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 25, fontface = "bold"),
        show_heatmap_legend=FALSE,
        use_raster = FALSE)
dev.off()

## plot legend ##
merge.df.dcast.temp = fread(paste0(indir,"/hyper_DHMR.5hmCG_mean_diff_old_minus_young_matrix.subclass_specific_top_1000.csv"),data.table=F,header=T)
rownames(merge.df.dcast.temp) = merge.df.dcast.temp[,1]
merge.df.dcast.temp  = merge.df.dcast.temp[,-1]
total_mC.diff_matrix = as.matrix(merge.df.dcast.temp)
inte.col = readRDS("../../input/01-youth/subclass_new.col_latest.rds")
inte.col = c(inte.col,c("CA2-FC-IG Glut"="#89C75F","L6b CTX Glut"="#0C727C","Lamp5 Lhx6 Gaba"="#90D5E4","PAL-STR Gaba-Chol" ="#00ae9d","DG-PIR Ex IMN"="#1d953f",
    "Vip Gaba"="#009ad6","L2/3 IT PIR-ENTl Glut"="#6E4B9E","L6 IT CTX Glut"="#AA0DFE","HPF CR Glut"='#e74c3c',"Pvalb chandelier Gaba"="#A6BDD7","STR D1 Sema5a Gaba" ="#B32851",
    "OB-mi Frmd7 Gaba"='#5AC2F1FF',"OB Trdn Gaba"="#e4c6d0","OB Meis2 Thsd7b Gaba"="#f9906f","VLMC NN"="#ffc773","Sst Chodl Gaba"="#88c4e8","STR Prox1 Lhx6 Gaba"="#eb7f54",
    "OT D3 Folh1 Gaba"="#815463","ABC NN" ="#253494","BAM NN" ="#FFFF00","Endo NN"="#d6ecf0","Peri NN"="#DEA0FD","Lymphoid NN"="#808080","SMC NN"="#bce672",'zeng'= 'lightgrey',
    "our"='lightgrey',"old"="lightgrey","young"="lightgrey","IT AON-TT-DP Glut"="#F6768E","LA-BLA-BMA-PA Glut"="#ff3300","COAa-PAA-MEA Barhl2 Glut"="#801dae"))
inte.col = inte.col[col_order]
annotation_df <- data.frame(
  subclass = factor(colnames(total_mC.diff_matrix),levels = col_order),
  stringsAsFactors = FALSE
)


ha1=HeatmapAnnotation(df = annotation_df, show_legend = TRUE,annotation_name_side = "right",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))

pdf("../../output/03-aging/03-DMRs_DHMRs/subclass_specific_top_DMRs_DHMRs_heatmap/heatmap_legend.max_0.2.pdf",width=8.2677,height=21.2598)
col_fun =colorRamp2(c(-0.2,-0.1,0,0.1,0.2), c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))
Heatmap(total_mC.diff_matrix,
        column_gap = unit(0, "points"),
        bottom_annotation = ha1,
        col = col_fun,
        cluster_rows=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names=F,
        show_row_names=F,
        row_names_side="left",
        column_names_side="bottom",
        column_names_rot = 75,
        column_dend_side = "top",
        border = TRUE,
        row_title = "Aging hyper-DHMRs (Aged > Young)",
        column_title = "5hmCG change",
        column_title_gp = gpar(fontsize = 25, fontface = "bold"),
        column_title_side = "top",
        column_title_rot = 0,
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 25, fontface = "bold"),
        heatmap_legend_param = list(title_position = "lefttop-rot",title = "diff (old - young)",legend_height = unit(4, "cm"),fontsize = 20),
        show_heatmap_legend=TRUE,
        use_raster = FALSE)
dev.off()