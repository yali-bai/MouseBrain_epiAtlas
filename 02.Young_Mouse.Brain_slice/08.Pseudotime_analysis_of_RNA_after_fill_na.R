library(scibetR)
library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(harmony)
library(monocle)

##### 01. pseudotime analysis #####
merged.seuratobj.sct <- readRDS("merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
gene.v <- rownames(merged.seuratobj.sct@assays$RNA$counts)
metainfo = merged.seuratobj.sct@meta.data

RNA_fill = readRDS("RNA_fill_na.rds")
RNA_fill.obj <- CreateSeuratObject(RNA_fill)
RNA_fill.obj <- NormalizeData(RNA_fill.obj, verbose = FALSE)
RNA_fill.obj <- ScaleData(RNA_fill.obj)
RNA_fill.obj <- FindVariableFeatures(RNA_fill.obj)
RNA_fill.obj <- RunPCA(RNA_fill.obj, verbose = FALSE)
RNA_fill.obj <- RunUMAP(RNA_fill.obj, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40)
RNA_fill.obj <- RunTSNE(RNA_fill.obj, dims=1:30, dim.embed=3)
RNA_fill.obj <- FindNeighbors(RNA_fill.obj, dims=1:30, reduction="pca")
RNA_fill.obj <- FindClusters(RNA_fill.obj)
length(which(rownames(RNA_fill.obj@meta.data) != rownames(metainfo)))
RNA_fill.obj$subclass <-  metainfo$subclass

RNA_fill.obj$select <- "N"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L2/3 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L4/5 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L5 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L6 IT CTX Glut")] <- "Y"
# decode comments if process CA Glut pseudotime
# RNA_fill.obj$select[which(RNA_fill.obj$subclass == "CA1-ProS Glut")] <- "Y"
# RNA_fill.obj$select[which(RNA_fill.obj$subclass == "CA2-FC-IG Glut")] <- "Y"
# RNA_fill.obj$select[which(RNA_fill.obj$subclass == "CA3 Glut")] <- "Y"

IT.obj <- subset(RNA_fill.obj, select == "Y")
IT.obj$layer <- IT.obj$subclass


## extract phenotypic information, that is, subclass information ##
expr_matrix <- as(as.matrix(IT.obj@assays$RNA$counts), 'sparseMatrix')
## extract phenotypic information to p_data(phenotype_data) ##
p_data <- IT.obj@meta.data[,-which(colnames(IT.obj@meta.data) == "subclass")]
## extract gene information ##
f_data <- data.frame(data.frame(gene_short_name = row.names(IT.obj@assays$RNA$counts), row.names = row.names(IT.obj@assays$RNA$counts)))
## The expression value matrix must: have the same number of columns as the phenoData has rows, and have the same number of rows as the featureData data frame has rows.
## create a new CellDataSet object (required) ##
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

monocle_cds <- newCellDataSet(expr_matrix,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

## estimate size factors and dispersions (required) ##
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
## Filtering low-quality cells ##
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1) 
print(head(fData(monocle_cds)))
expressed_genes <- row.names(subset(fData(monocle_cds), num_cells_expressed >= 10)) 

## choose genes that define a cell's progress ##
Idents(IT.obj) <- 'layer'
diff <- FindAllMarkers(IT.obj,assay = "RNA", slot = "counts",logfc.threshold = 0.5,only.pos = T,min.pct = 0.2)
diff <- subset(diff,!grepl(pattern = 'RP[LS]',gene))
diff <- subset(diff,!grepl(pattern = 'MT-',gene))
head(diff)


ordergene <- row.names(diff) 
monocle_cds <- setOrderingFilter(monocle_cds, ordergene)
#plot_ordering_genes(monocle_cds)

## reduce data dimensionality ##
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                        method = "DDRTree")
## order cells along the trajectory ##
monocle_cds <- orderCells(monocle_cds)

ordergene <- intersect(monocle_cds@featureData@data$gene_short_name, ordergene)
## order cells along the trajectory ##
monocle_cds <- orderCells(monocle_cds)

## perform differential gene expression test as a way to extract the genes that distinguish them ##
sig_diff.genes <- subset(diff,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))
diff_test <- differentialGeneTest(monocle_cds[sig_diff.genes,], cores = 1, 
                              fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))

## calculate subclass mean of candidate subclasses ##
subclass.df <- data.frame(matrix(NA,nrow=1122,ncol=1))
for(cl in unique(IT.obj$layer)){
    temp.obj <- subset(IT.obj, layer == cl)
    subclass.df <- cbind(subclass.df,rowMeans(temp.obj@assays$RNA$counts))
}
subclass.df <- subclass.df[,-1]
colnames(subclass.df) <- unique(IT.obj$layer)

Pseudotime_gene.df <- subclass.df[rownames(diff_test[which(diff_test$use_for_ordering),]),]
Pseudotime_gene.df = Pseudotime_gene.df[,c("L2/3 IT CTX Glut","L4/5 IT CTX Glut","L5 IT CTX Glut","L6 IT CTX Glut")]
#saveRDS(Pseudotime_gene.df,file="Pseudotime_genes_RNA_expr.IT_Glut.rds")

IT_sort_decreasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] < Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] < Pseudotime_gene.df[,3] & Pseudotime_gene.df[,3] < Pseudotime_gene.df[,4],]                   
IT_sort_increasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] > Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] > Pseudotime_gene.df[,3] & Pseudotime_gene.df[,3] > Pseudotime_gene.df[,4],]
save(IT_sort_decreasing.df,IT_sort_increasing.df,file="IT_pseudotime_genes.RData")
# CA_sort_decreasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] < Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] < Pseudotime_gene.df[,3],]
# CA_sort_increasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] > Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] > Pseudotime_gene.df[,3],]
# order_Pseudotime_gene.df <- Pseudotime_gene.df[order(Pseudotime_gene.df [,1]),]
# save(CA_sort_increasing.df,CA_sort_decreasing.df,file="CA_pseudotime_genes.RData")

##### plot candidate pseudotime genes #####
## packages
library(Seurat)
now_lib <- .libPaths()
.libPaths(c(now_lib,"/share/home/zhangac/anaconda3/envs/Seurat/lib/R/library"))
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(MuDataSeurat)
library(grid)
library(cowplot)
library(reshape2)


order.v = readRDS("../03.data/04.config_files/subclass_order.rds")

paired_sampleinfo = read.csv("../03.data/02.metainfo/02.Young_Mouse.Brain_slice/RNA_DNA_match_name_QC_class_label_young.brain_slice.add_celltype.csv",header=T)


##### spatial color #####
loc.df <- read.csv("../03.data/03.download_data/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
# Zhuang MERFISH metainfo
rownames(loc.df) <- loc.df$cell_label
TSNE_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'tsne')[Cells(merged.seuratobj.sct),])  # extract final tsne loci
spatial_zhuang_merfish_loci.df <- as.data.frame(TSNE_RNA[which(is.na(str_match(rownames(TSNE_RNA),"TSO"))),])
spatial_zhuang_merfish_loci.df$spatial_1 <- loc.df[rownames(spatial_zhuang_merfish_loci.df),"x"]
spatial_zhuang_merfish_loci.df$spatial_2 <- loc.df[rownames(spatial_zhuang_merfish_loci.df),"y"]

##### metainfo
metainfo = merged.seuratobj.sct@meta.data
TSNE_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  # extract final tsne loci

load("DNA_fill_na.filter_by_total_QC.RData")
RNA_fill = readRDS("RNA_fill_na.rds")

##### set color
colors <- colorRampPalette(c("#253494", "#225EA8", "#1D91C0", "#41B6C4","#FFED6F","#FFFF33","#FFFF00")) # 从蓝色到黄色V2
color_gradient <- colors(75)

colors <- colorRampPalette(c("#5614b0", "#7438bd","#915cca","#af80d6","#cca4e3","#F4C800","#FFED6F","#FFFF00")) ## purple to yellow
color_gradient_RNA <- colors(75)

marker_plot <- function(cell_name,methyl_matrix,name,gene){
        if((name == "5hmC_CGN") || (name == "5hmC_CHN") || (name == "5hmCG") || (name == "5hmCH")){
            str_re="joint5hmC"
            rownames(paired_sampleinfo) = paired_sampleinfo$hmC
            coln = "hmC"
        }else{
            str_re="joint5mC"
            rownames(paired_sampleinfo) = paired_sampleinfo$mC
            coln = "mC"
        }
        plot_df <- as.data.frame(cbind(as.numeric(methyl_matrix[gene,]),as.numeric(TSNE_RNA[paired_sampleinfo[colnames(methyl_matrix),"RNA"],"umap_1"]),as.numeric(TSNE_RNA[paired_sampleinfo[colnames(methyl_matrix),"RNA"],"umap_2"])))
        colnames(plot_df) <- c("methyl","umap_1","umap_2")
        rownames(plot_df) <- colnames(methyl_matrix)
        #plot_df <- plot_df[str_replace(cell_name,"RNA",str_re),]
        rownames(paired_sampleinfo) = paired_sampleinfo$RNA
        plot_df <- plot_df[paired_sampleinfo[cell_name,coln],]
        plot_df <- na.omit(plot_df)
        bigger_than_zero.idx <- which(as.numeric(plot_df$methyl) != 0)
        equal_to_zero.idx <- which(as.numeric(plot_df$methyl) == 0)
        bigger_than_zero.df <- plot_df[bigger_than_zero.idx,]
        equal_to_zero.df <- plot_df[equal_to_zero.idx,]
        ggplot() + 
            geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") + 
            geom_point(data = equal_to_zero.df, mapping = aes(umap_1, umap_2, color=methyl), size = 6) +
            geom_point(data = bigger_than_zero.df, mapping = aes(umap_1, umap_2, color=methyl), size = 6) +
            scale_color_gradientn(colors = color_gradient,breaks=c(0,0.5,1))+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) +
            labs(title = name, x = "spatial_1", y = "spatial_2") 
}

marker_plot_limits <- function(cell_name,methyl_matrix,name,gene,min,max){
        if((name == "5hmC_CGN") || (name == "5hmC_CHN") || (name == "5hmCG") || (name == "5hmCH")){
            str_re="joint5hmC"
            rownames(paired_sampleinfo) = paired_sampleinfo$hmC_SampleID
            coln = "hmC_SampleID"
        }else{
            str_re="joint5mC"
            rownames(paired_sampleinfo) = paired_sampleinfo$mC_SampleID
            coln = "mC_SampleID"
        }
        plot_df <- as.data.frame(cbind(as.numeric(methyl_matrix[gene,]),as.numeric(TSNE_RNA[paired_sampleinfo[colnames(methyl_matrix),"RNA_SampleID"],"umap_1"]),as.numeric(TSNE_RNA[paired_sampleinfo[colnames(methyl_matrix),"RNA_SampleID"],"umap_2"])))
        colnames(plot_df) <- c("methyl","umap_1","umap_2")
        rownames(plot_df) <- colnames(methyl_matrix)
        rownames(paired_sampleinfo) = paired_sampleinfo$RNA_SampleID
        plot_df <- plot_df[paired_sampleinfo[cell_name,coln],]
        plot_df <- na.omit(plot_df)
        bigger_than_zero.idx <- which(as.numeric(plot_df$methyl) != 0)
        equal_to_zero.idx <- which(as.numeric(plot_df$methyl) == 0) 
        plot_df$methyl[which(plot_df$methyl < min)] = min
        bigger_than_zero.df <- plot_df[bigger_than_zero.idx,]
        equal_to_zero.df <- plot_df[equal_to_zero.idx,]
        ggplot() + 
            geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") + 
            geom_point(data = equal_to_zero.df, mapping = aes(umap_1, umap_2, color=methyl), size = 6) +
            geom_point(data = bigger_than_zero.df, mapping = aes(umap_1, umap_2, color=methyl), size = 6) +
            scale_color_gradientn(  
                limits = c(min, max), 
                breaks = c(min,max),
                colors = color_gradient,
                na.value = color_gradient[length(color_gradient)])+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30),
                legend.title = element_blank()) +
            labs(title = name, x = "spatial_1", y = "spatial_2") 
}

spatial_loci <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  
spatial_loci = as.matrix(spatial_loci*-1)
merged.seuratobj.sct@reductions$umap@cell.embeddings <- spatial_loci
TSNE_RNA = TSNE_RNA*-1
spatial_zhuang_merfish_loci.df = spatial_zhuang_merfish_loci.df*-1
spatial_zhuang_merfish_loci.df$umap_1 = spatial_zhuang_merfish_loci.df$spatial_1
spatial_zhuang_merfish_loci.df$umap_2 = spatial_zhuang_merfish_loci.df$spatial_2


limits_marker_RNA_cluster_plot_neuron_type <- function(marker_type, type, feature, cl, gene, pos_neg,max_vector){
                if(cl == "L2/3 IT CTX Glut"){
                    inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = 'lightgrey', "Astro-TE NN" = 'lightgrey',
                        "L2/3 IT CTX Glut" = "#ed1941", "Sst Gaba" = 'lightgrey', "L6 IT CTX Glut" = "#89288F", 
                        "OPC NN" = 'lightgrey', "L5 ET CTX Glut" = 'lightgrey', "Microglia NN" = 'lightgrey',  
                        "Peri NN" = 'lightgrey', "CA2-FC-IG Glut" = 'lightgrey', "Lamp5 Gaba" = 'lightgrey', 
                        "Endo NN" = 'lightgrey', "DG Glut" = 'lightgrey', "L4/5 IT CTX Glut" = "#009ad6", 
                        "Vip Gaba" = 'lightgrey', "L6 CT CTX Glut" = 'lightgrey', "Pvalb Gaba" = 'lightgrey', 
                        "CLA-EPd-CTX Car3 Glut" = 'lightgrey', "L2/3 IT RSP Glut" = 'lightgrey', "CA3 Glut" = 'lightgrey',
                        "Ependymal NN"='lightgrey',"CHOR NN"='lightgrey',"HPF CR Glut"='lightgrey',"CA1-ProS Glut"='lightgrey',
                        "L5 NP CTX Glut"='lightgrey',"SMC NN"='lightgrey',"VLMC NN"='lightgrey',"L4 RSP-ACA Glut"='lightgrey',"DG-PIR Ex IMN"='lightgrey',
                        "L5 IT CTX Glut" = "#FFB300", "L6b CTX Glut" = 'lightgrey',"Lamp5 Lhx6 Gaba"='lightgrey',"RHP-COA Ndnf Gaba"='lightgrey',
                        "Sncg Gaba"='lightgrey',"Sst Chodl Gaba"='lightgrey',"STR D2 Gaba"='lightgrey',"OB-STR-CTX Inh IMN"="lightgrey")
                    Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass
                    rm.idx = c(which(levels(merged.seuratobj.sct) == "L2/3 IT CTX Glut"),
                               which(levels(merged.seuratobj.sct) == "L6 IT CTX Glut"),
                               which(levels(merged.seuratobj.sct) == "L4/5 IT CTX Glut"),
                               which(levels(merged.seuratobj.sct) == "L5 IT CTX Glut"))
                    all_subclass_plot <- DimPlot(merged.seuratobj.sct, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =6,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("L2/3 IT CTX Glut","L4/5 IT CTX Glut","L5 IT CTX Glut","L6 IT CTX Glut",levels(merged.seuratobj.sct)[-rm.idx]))+
                        theme(  
                            axis.title.x = element_blank(), 
                            axis.text.x = element_blank(),  
                            axis.ticks.x = element_blank(), 
                            axis.title.y = element_blank(), 
                            axis.text.y = element_blank(),  
                            axis.ticks.y = element_blank(), 
                            line = element_blank()
                        )
                    cell_n.v <- rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "L2/3 IT CTX Glut"),])
                    cell_n.v <-c(cell_n.v,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "L4/5 IT CTX Glut"),]))
                    cell_n.v <-c(cell_n.v,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "L5 IT CTX Glut"),]))
                    cell_n.v <-c(cell_n.v,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "L6 IT CTX Glut"),]))
                }else{
                    inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = 'lightgrey', "Astro-TE NN" = 'lightgrey',
                        "L2/3 IT CTX Glut" = 'lightgrey', "Sst Gaba" = 'lightgrey', "L6 IT CTX Glut" = 'lightgrey', 
                        "OPC NN" = 'lightgrey', "L5 ET CTX Glut" = 'lightgrey', "Microglia NN" = 'lightgrey',  
                        "Peri NN" = 'lightgrey', "CA2-FC-IG Glut" = "#f173ac", "Lamp5 Gaba" = 'lightgrey', 
                        "Endo NN" = 'lightgrey', "DG Glut" = 'lightgrey', "L4/5 IT CTX Glut" = 'lightgrey', 
                        "Vip Gaba" = 'lightgrey', "L6 CT CTX Glut" = 'lightgrey', "Pvalb Gaba" = 'lightgrey', 
                        "CLA-EPd-CTX Car3 Glut" = 'lightgrey', "L2/3 IT RSP Glut" = 'lightgrey', "CA3 Glut" = "#00538A",
                        "Ependymal NN"='lightgrey',"CHOR NN"='lightgrey',"HPF CR Glut"='lightgrey',"CA1-ProS Glut"="#2585a6",
                        "L5 NP CTX Glut"='lightgrey',"SMC NN"='lightgrey',"VLMC NN"='lightgrey',"L4 RSP-ACA Glut"='lightgrey',"DG-PIR Ex IMN"='lightgrey',
                        "L5 IT CTX Glut" = 'lightgrey', "L6b CTX Glut" = 'lightgrey',"Lamp5 Lhx6 Gaba"='lightgrey',"RHP-COA Ndnf Gaba"='lightgrey',
                        "Sncg Gaba"='lightgrey',"Sst Chodl Gaba"='lightgrey',"STR D2 Gaba"='lightgrey',"OB-STR-CTX Inh IMN"="lightgrey")
                        Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass
                        rm.idx = c(which(levels(merged.seuratobj.sct) == "CA3 Glut"),
                               which(levels(merged.seuratobj.sct) == "CA1-ProS Glut"),
                               which(levels(merged.seuratobj.sct) == "CA2-FC-IG Glut"))
                        all_subclass_plot <- DimPlot(merged.seuratobj.sct, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =6,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("CA1-ProS Glut","CA2-FC-IG Glut","CA3 Glut",levels(merged.seuratobj.sct)[-rm.idx]))+
                            theme(  
                                axis.title.x = element_blank(), 
                                axis.text.x = element_blank(),  
                                axis.ticks.x = element_blank(), 
                                axis.title.y = element_blank(), 
                                axis.text.y = element_blank(),  
                                axis.ticks.y = element_blank(), 
                                line = element_blank()
                            )
                        cell_n.v <- rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "CA1-ProS Glut"),])
                        cell_n.v <-c(cell_n.v,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "CA2-FC-IG Glut"),]))
                        cell_n.v <-c(cell_n.v,rownames(merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$subclass == "CA3 Glut"),]))
                }
                plot_df <- as.data.frame(cbind(as.numeric(RNA_fill[strsplit(gene,"\\.")[[1]][1],which(!is.na(str_match(colnames(RNA_fill),"TSO")))]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))],"umap_1"]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))],"umap_2"])))
                colnames(plot_df) <- c("RNA","umap_1","umap_2")
                rownames(plot_df) <- colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))]

                plot_df <- plot_df[cell_n.v,]
                plot_df <- na.omit(plot_df)
                plot_df[which(plot_df$RNA < max_vector[1]),"RNA"] = max_vector[1]
                spatial_zhuang_merfish_loci.df <- TSNE_RNA[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO"))),]

                RNA_plot <- ggplot() + 
                geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") +
                geom_point(data = plot_df, mapping = aes(umap_1, umap_2, color=RNA), size = 6) +
                scale_color_gradientn(  
                    limits = c(max_vector[1], max_vector[2]), 
                    breaks = seq(max_vector[1], max_vector[2],1),
                    colors = color_gradient_RNA,
                    na.value = color_gradient_RNA[length(color_gradient_RNA)])+
                theme_minimal() +
                theme(text = element_text(face="bold", size=30),panel.grid = element_blank(),axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_text(face="bold", size=30),
                    axis.title.y = element_text(face="bold", size=30))+
                labs(title = "RNA", x = "spatial_1", y = "spatial_2")

                zhuang_plot_df <- as.data.frame(cbind(as.numeric(RNA_fill[strsplit(gene,"\\.")[[1]][1],which(is.na(str_match(colnames(RNA_fill),"TSO")))]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))],"umap_1"]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))],"umap_2"])))
                colnames(zhuang_plot_df) <- c("RNA","umap_1","umap_2")
                rownames(zhuang_plot_df) <- colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"TSO")))]
                zhuang_plot_df <- zhuang_plot_df[cell_n.v,]
                zhuang_plot_df <- na.omit(zhuang_plot_df)
                zhuang_plot_df[which(zhuang_plot_df$RNA < max_vector[3]),"RNA"] = max_vector[3]
                zhuang_RNA_plot <- ggplot() + 
                geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") +
                geom_point(data = zhuang_plot_df, mapping = aes(umap_1, umap_2, color=RNA), size = 6) +
                scale_color_gradientn(  
                    limits = c(max_vector[3], max_vector[4]), 
                    breaks = seq(max_vector[3], max_vector[4],1),
                    colors = color_gradient_RNA,
                    na.value = color_gradient_RNA[length(color_gradient_RNA)])+
                theme_minimal() +
                theme(text = element_text(face="bold", size=30),panel.grid = element_blank(),axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_text(face="bold", size=30),
                    axis.title.y = element_text(face="bold", size=30))+
                labs(title = "RNA.Zhuang", x = "spatial_1", y = "spatial_2") 

                if(feature == "genebody"){
                    ## mCG
                    mC_CGN = gene_mCG_fill
                                    
                    ## hmCG
                    hmC_CGN = gene_hmCG_fill
                    
                    ## mCG_hmCG
                    hmCG_mCG = gene_hmCG_mCG_fill                                 
                }else if(feature == "promoter"){
                    ## mCG
                    mC_CGN = promoter_mC_CGN_fill
                
                    ## hmCG
                    hmC_CGN = promoter_hmC_CGN_fill

                    ## mCG_hmCG
                    hmCG_mCG = promoter_hmCG_mCG_fill
                }

                for(i in 1:dim(mC_CGN)[2]){
                    mC_CGN[which(mC_CGN[,i] < 0),i] <- 0
                }


                mCG_plot <- marker_plot_limits(cell_n.v,mC_CGN,"5mCG",strsplit(gene,"\\.")[[1]][1],max_vector[5],max_vector[6])
                hmCG_plot <- marker_plot_limits(cell_n.v,hmC_CGN,"5hmCG",strsplit(gene,"\\.")[[1]][1],max_vector[7],max_vector[8])
                mCG_hmCG_plot <- marker_plot_limits(cell_n.v,hmCG_mCG,"5mCG+5hmCG",strsplit(gene,"\\.")[[1]][1],max_vector[9],max_vector[10])

                plot_grid(  
                    plotlist = list(  
                        plot_grid(all_subclass_plot, RNA_plot, zhuang_RNA_plot, ncol = 3, rel_widths = c(8,8,8)),  
                        plot_grid(mCG_plot, hmCG_plot, mCG_hmCG_plot, ncol = 3, rel_widths = c(8,8,8))
                    ),  
                    nrow = 2
                )

}


plot <- limits_marker_RNA_cluster_plot_neuron_type("RNA","subclass","genebody",'L2/3 IT CTX Glut','ENSMUSG00000001985.9',"positive",c(0,2,0,4,0.35,0.65,0.25,0.4,0.5,0.9))    
pdf("L2_3_IT_CTX_Glut.ENSMUSG00000001985.pdf",width = 50,height = 25)
print(plot)
dev.off()

plot <- limits_marker_RNA_cluster_plot_neuron_type("RNA","subclass","genebody",'CA3 Glut','ENSMUSG00000032537.15',"positive",c(0,2,0,4,0.35,0.65,0.25,0.4,0.5,0.9))    
pdf("CA3_Glut.ENSMUSG00000032537.pdf",width = 50,height = 25)
print(plot)
dev.off()



