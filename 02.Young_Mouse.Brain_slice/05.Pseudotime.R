##### 01.import packages #####
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
## https://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories ##


##### 02.change working path #####
# setwd("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/MERFISH/20240902")

##### 03.read subclass information #####
merged.seuratobj.sct <- readRDS("../output/02-slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
gene.v <- rownames(merged.seuratobj.sct@assays$RNA$counts)
metainfo = merged.seuratobj.sct@meta.data

##### 04.monocle analysis process #####
load("../output/02-slice/RNA_DNA_fill_na.20240925.RData")
## pretreatment ##
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
RNA_fill.obj$subclass <-  metainfo$subclass ## assign subclass information

## select candidate cell type ##
RNA_fill.obj$select <- "N"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L2/3 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L4/5 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L5 IT CTX Glut")] <- "Y"
RNA_fill.obj$select[which(RNA_fill.obj$subclass == "L6 IT CTX Glut")] <- "Y"
# run the following lines if CA Glut was analysised
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
                              expressionFamily = negbinomial.size()) # raw read counts

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

#ordergene <- intersect(monocle_cds@featureData@data$gene_short_name, ordergene)
#monocle_cds <- orderCells(monocle_cds)

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
IT_sort_decreasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] < Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] < Pseudotime_gene.df[,3] & Pseudotime_gene.df[,3] < Pseudotime_gene.df[,4],]                   
IT_sort_increasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] > Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] > Pseudotime_gene.df[,3] & Pseudotime_gene.df[,3] > Pseudotime_gene.df[,4],]
save(IT_sort_decreasing.df,IT_sort_increasing.df,file="../output/02-slice/pseudotime/IT_pseudotime_genes.RData")

# run the following lines if CA Glut was analysised
# CA_sort_decreasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] < Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] < Pseudotime_gene.df[,3],]
# CA_sort_increasing.df <- Pseudotime_gene.df[Pseudotime_gene.df[,1] > Pseudotime_gene.df[,2] & Pseudotime_gene.df[,2] > Pseudotime_gene.df[,3],]
# order_Pseudotime_gene.df <- Pseudotime_gene.df[order(Pseudotime_gene.df [,1]),]
# save(CA_sort_increasing.df,CA_sort_decreasing.df,file="../output/02-slice/pseudotime/CA_pseudotime_genes.RData")
