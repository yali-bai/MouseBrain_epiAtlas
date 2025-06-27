##### 01. import packages #####
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

##### 03. set subclass color #####
inte.col = readRDS("../../04.data/04.config_files/color.subclass.rds")

##### 04. subclass order #####
subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")

##### 05. decide DHMRs region order and intersected gene order #####
##### all hyper DHMRs intersected genes ##### 
## filter top1000 hyper DHMRs later ##
subclass = c()
for(dir in c("all")){
    result = data.frame()
    for(cl in subclass_order){
        subclass_filename = str_replace_all(str_replace_all(cl,"/","_")," ","_")
        file_exists <- file.exists(paste0("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/08.DHMRs_related_gene_GO_analysis/20250217_update.filter_before_calculating_p_value_adjust/subclass/",dir,"/",subclass_filename,".",dir,"_hyper_DHMRs.bed.intersection.20250410.wo.bed"))
        if (file_exists) {
            gene_df = read.table(paste0("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/08.DHMRs_related_gene_GO_analysis/20250217_update.filter_before_calculating_p_value_adjust/subclass/",dir,"/",subclass_filename,".",dir,"_hyper_DHMRs.bed.intersection.20250410.wo.bed"),header=F,sep="\t")
            result = rbind(result,gene_df)
            subclass = c(subclass,rep(cl,dim(gene_df)[1]))
        }
    }
}
result$subclass = subclass
dim(result)

##### retain high express and high CV gene, but do not restrict the intersection length of DHMRs and gene #####
colnames(result) = c("DHMR_chr","DHMR_start","DHMR_end","gene_chr","gene_start","gene_end","geneid","length","subclass")
## open if you want filter intersection length ##
#result = result[result$length > 500,]
result$DHMR_region = paste0(result$DHMR_chr,"_",result$DHMR_start,"_",result$DHMR_end)
result$geneid = unlist(lapply(result$geneid, function(x) strsplit(x,"\\.")[[1]][1]))
## filter high expressed and high CV genes ##
filter_gene = read.table("../../04.data/05.intermediate_files/01.RNA/02.Aging_Mouse/HighExprHighCV_genes.txt",header=F)
result = result[result$geneid %in% filter_gene$V1,]
head(result)

##### significant DHMR #####
DHMR_result=readRDS(paste0(outdir,"/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds"))
head(DHMR_result)

##### top 1000 hyper #####
result$uniq = paste0(result$subclass,";",result$DHMR_region)
DHMR_result$uniq = paste0(DHMR_result$cluster,";",DHMR_result$chrom)
result$diff = DHMR_result[match(result$uniq,DHMR_result$uniq),"diff"]
result$subclass = factor(result$subclass,levels = subclass_order)
head(result)
result.sorted_by_diff_in_each_subclass = result %>%
    filter(diff > 0) %>%
    group_by(subclass) %>%
    arrange(desc(diff)) %>%
    slice_head(n=1000)
head(result.sorted_by_diff_in_each_subclass)

dim(result.sorted_by_diff_in_each_subclass)

## sort by subclass_order and diff ##
result.sorted_by_diff_in_each_subclass.sort = result.sorted_by_diff_in_each_subclass[order(result.sorted_by_diff_in_each_subclass$subclass,-result.sorted_by_diff_in_each_subclass$diff),]
head(result.sorted_by_diff_in_each_subclass.sort)

##### plot order #####
DHMR_region_order = rev(result.sorted_by_diff_in_each_subclass.sort$DHMR_region)
gene_order = rev(result.sorted_by_diff_in_each_subclass.sort$geneid)


##### 06. plot DHMR heatmap #####
## read DHMR region methylation diff (old - young) data ##
merge.df.dcast = readRDS(paste0(indir,"merge.df.dcast.rds"))

total_mC.dcast = merge.df.dcast[which(merge.df.dcast$type == "5hmCG"),]
total_mC.dcast_subset = total_mC.dcast[,c("segment","subclass","diff")]
total_mC.diff_matrix = dcast(total_mC.dcast_subset,segment~subclass)

rownames(total_mC.diff_matrix) = total_mC.diff_matrix$segment
total_mC.diff_matrix = total_mC.diff_matrix[,-1]
total_mC.diff_matrix = total_mC.diff_matrix[DHMR_region_order,]
total_mC.diff_matrix = total_mC.diff_matrix[,intersect(subclass_order,colnames(total_mC.diff_matrix))]

total_mC.diff_matrix[1:5,1:5]

total_mC.diff_matrix = as.matrix(total_mC.diff_matrix)
length(inte.col)
inte.col = inte.col[colnames(total_mC.diff_matrix)]
length(inte.col)

ha1=rowAnnotation(subclass = factor(rev(result.sorted_by_diff_in_each_subclass.sort$subclass),levels=subclass_order),show_legend = FALSE,annotation_name_side = "top",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))
ha2=HeatmapAnnotation(subclass = colnames(total_mC.diff_matrix),show_legend = TRUE,annotation_name_side = "right",col = list(subclass=inte.col),annotation_legend_param = list(direction = "horizontal", title_position = "topleft",fontsize = 20))

pdf("../../output/03.Aging_Mouse/03-DMRs_DHMRs/top1000_hyper_DHMR_intersected_with_highexpressed_highCV_gene.DHMRs_5hmCG_mean_diff_old_minus_young_matrix.sorted_by_subclass_and_diff.heatmap.20250612.pdf",width=7,height=10.2518)
col_fun =colorRamp2(c(-0.15,-0.075,0,0.075,0.15), c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3")) 
Heatmap(total_mC.diff_matrix,
        column_gap = unit(0, "points"),
        top_annotation = ha2,
        left_annotation = ha1,
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
        use_raster = FALSE)
dev.off()

##### intersect gene with top1 subclass marker labeled #####
gene_expression = read.csv(paste0(indir,"/Joint_Cabernet_RNA_data_old_young_diff_total_result_corrected_all_gene.csv",header=T))
colnames(gene_expression) = c("geneid","level","subclass","young_mean","old_mean","diff","young_0_ratio","old_0_ratio")
head(gene_expression)

result_expression = gene_expression[gene_expression$geneid %in% gene_order,]
head(result_expression)
result_expression$subclass = factor(result_expression$subclass,levels=subclass_order)
result_expression$log2FoldChange = log(result_expression$old_mean/result_expression$young_mean)
head(result_expression)

mm10_gene_metainfo = read.table("../../04.data/01.ref/mm10.genes_duplicated.bed",header=T,sep="\t")
mm10_gene_metainfo$gene_id = unlist(lapply(mm10_gene_metainfo$gene_id, function(x) strsplit(x,"\\.")[[1]][1]))
rownames(mm10_gene_metainfo) = mm10_gene_metainfo$gene_id
head(mm10_gene_metainfo)


result_expression$genename = mm10_gene_metainfo[result_expression$geneid,"gene_name"]
head(result_expression)

result_expression.temp = result_expression[,c("subclass","geneid","young_mean")]
head(result_expression.temp)
result_expression.young = dcast(result_expression.temp,geneid~subclass, value.var = "young_mean")
rownames(result_expression.young) = result_expression.young$geneid
result_expression.young = result_expression.young[,-1]
head(result_expression.young)

result_expression.temp = result_expression[,c("subclass","geneid","old_mean")]
result_expression.aged = dcast(result_expression.temp,geneid~subclass)
rownames(result_expression.aged) = result_expression.aged$geneid
result_expression.aged = result_expression.aged[,-1]
head(result_expression.aged)

meanscore <- function(x) {
  return((x / max(x, na.rm = TRUE)))
}

top1_df = read.csv("../../04.data/05.intermediate_files/01.RNA/02.Aging_Mouse/top1000_hyper_DHMRs_intersected_gene.intersected_subclass_marker.top1_log2FC_in_subclass_markers.csv",header=T)
top1_df$loci = paste0(top1_df$subclass,";",top1_df$geneid)
#head(top1_df[top1_df$subclass == "Astro-TE NN",])

gene_order_df = result.sorted_by_diff_in_each_subclass.sort[which(result.sorted_by_diff_in_each_subclass.sort$DHMR_region %in% DHMR_region_order),]
gene_order_df$loci = paste0(gene_order_df$subclass,";",gene_order_df$geneid)
#head(gene_order_df[gene_order_df$subclass == "Astro-TE NN" & gene_order_df$geneid == "ENSMUSG00000024411",])

dim(result.sorted_by_diff_in_each_subclass.sort)
dim(gene_order_df)

index = match(top1_df$loci,rev(gene_order_df$loci))

top1_markers = top1_df$genename

result.zscore = apply(as.matrix(result_expression.young), 1, meanscore)
result.zscore[1:5,1:5]
result.zscore = result.zscore[,gene_order]
result.zscore[1:5,1:5]
colnames(result.zscore) = mm10_gene_metainfo[colnames(result.zscore),"gene_name"]
result.zscore[1:5,1:5]

col_fun =colorRamp2(c(0,0.55,0.7,1), c("#d3d3d3","#e9e9e9","#A400D3","#5A0099")) #v20
ann_colors=list(subclass=inte.col[rownames(result.zscore)])
annotation_col = data.frame(subclass=rownames(result.zscore))
annotation_col$subclass = factor(annotation_col$subclass,levels=subclass_order)

phtm1 <- ComplexHeatmap::pheatmap(t(result.zscore), use_raster=F, border_color = NA, fontsize = 5, color = colorRamp2(c(0,0.55,0.7,1), c("#d3d3d3","#e9e9e9","#A400D3","#5A0099")), 
                                 labels_row = NULL, cellwidth = 10, cluster_cols = F, cluster_rows = F,show_colnames = FALSE,
                                 annotation_colors=ann_colors, annotation_col=annotation_col)

ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=index, 
                                                                  labels=top1_markers))
markerGene_heatmap_zeng <- phtm1 + ha
options(repr.plot.width = 15, repr.plot.height = 9, repr.plot.res = 100)
markerGene_heatmap_zeng
pdf("top1000_hyper_DHMR_intersected_high_express_and_highCV_gene.RNA_expression.max_score_by_gene.color_v20.label_top1_subclass_marker.young.250612.pdf",width=9.1212,height=10.2518) #21.2598
print(markerGene_heatmap_zeng)
dev.off()

top1_markers = top1_df$genename

result.zscore = apply(as.matrix(result_expression.aged), 1, meanscore)
result.zscore[1:5,1:5]
result.zscore = result.zscore[,gene_order]
result.zscore[1:5,1:5]
colnames(result.zscore) = mm10_gene_metainfo[colnames(result.zscore),"gene_name"]
result.zscore[1:5,1:5]

col_fun =colorRamp2(c(0,0.55,0.7,1), c("#d3d3d3","#e9e9e9","#A400D3","#5A0099")) #v20
ann_colors=list(subclass=inte.col[rownames(result.zscore)])
annotation_col = data.frame(subclass=rownames(result.zscore))
annotation_col$subclass = factor(annotation_col$subclass,levels=subclass_order)

phtm1 <- ComplexHeatmap::pheatmap(t(result.zscore), use_raster=F, border_color = NA, fontsize = 5, color = colorRamp2(c(0,0.55,0.7,1), c("#d3d3d3","#e9e9e9","#A400D3","#5A0099")), 
                                 labels_row = NULL, cellwidth = 10, cluster_cols = F, cluster_rows = F,show_colnames = FALSE,
                                 annotation_colors=ann_colors, annotation_col=annotation_col)

ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=index, 
                                                                  labels=top1_markers))
markerGene_heatmap_zeng <- phtm1 + ha
options(repr.plot.width = 15, repr.plot.height = 9, repr.plot.res = 100)
markerGene_heatmap_zeng
pdf("top1000_hyper_DHMR_intersected_high_express_and_highCV_gene.RNA_expression.max_score_by_gene.color_v20.label_top1_subclass_marker.aged.250612.pdf",width=9.1212,height=10.2518) #21.2598
print(markerGene_heatmap_zeng)
dev.off()


