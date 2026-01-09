######## Integration ########
library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
library(stringr)
options(future.globals.maxSize = 8000 * 1024^2)

meta <- read.delim('01-meta.txt',head=F,stringsAsFactors=F)
topdim <- 1:30 

seuratObj.list <- lapply(meta$V1, function(x){
  readRDS(sprintf("rds/%s.rds", x))
})

Joint_Cabernet_rds = readRDS("input/Joint_Cabernet.rds")
non_zero_ratio <- rowMeans(Joint_Cabernet_rds > 0)
genes = names(non_zero_ratio[non_zero_ratio <= 0.8])

mm10_geneID <- read.table("../../../03.data/01.ref/mm10.genes.bed",header=F)
colnames(mm10_geneID) = c("chr","start","end","geneid","gene_name","gene_type")
mm10_geneID = mm10_geneID[-which(mm10_geneID$chr %in% c("chrX","chrY","chrM")),]
mm10_geneID = mm10_geneID[mm10_geneID$gene_type == "protein_coding",]
mm10_geneID$geneid_without_version = unlist(lapply(mm10_geneID$geneid, function(x) strsplit(x,'\\.')[[1]][1]))

genes = intersect(genes,mm10_geneID$geneid)
common_genes <- Reduce(intersect, lapply(seuratObj.list, function(obj) {
  rownames(GetAssayData(obj, assay = "RNA",layer="counts"))  
}))
common_genes = intersect(common_genes,genes)

seuratObj.list <- lapply(seuratObj.list, function(obj) {
  subset(obj, features = common_genes)
})

seuratObj.list <- lapply(X = seuratObj.list, FUN = SCTransform, method = "glmGamPoi")

top.features <- SelectIntegrationFeatures(object.list = seuratObj.list, nfeatures = 3000)

seuratObj.list <- PrepSCTIntegration(object.list = seuratObj.list, anchor.features = top.features)
seuratObj.list <- lapply(X = seuratObj.list, FUN = RunPCA, features = top.features)

integrated.anchors <- FindIntegrationAnchors(object.list = seuratObj.list, normalization.method = "SCT", 
    anchor.features = top.features, dims = topdim, reduction = "rpca", k.anchor = 10) # k.anchor = 20
integrated <- IntegrateData(anchorset = integrated.anchors, normalization.method = "SCT", dims = topdim)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
# integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunTSNE(integrated, dims=topdim, dim.embed=3) 
integrated <- RunUMAP(integrated, dims=topdim, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40) 
integrated <- FindNeighbors(integrated, dims=topdim, reduction="pca") 
integrated <- FindClusters(integrated, resolution=1.8) 
integrated@misc$geneName <- seuratObj.list[[1]]@misc$geneName

# Normalize RNA data for visualization purposes
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)

integrated <- JoinLayers(integrated)
DefaultAssay(integrated) <- "integrated"
saveRDS(integrated, "integration_4.rds")
