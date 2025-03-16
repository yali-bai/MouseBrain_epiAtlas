#########    All "our" in the following code refers to Joint Cabernet.

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

######## Integration ########
library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
options(future.globals.maxSize = 8000 * 1024^2)
# plan(strategy = "multicore", workers = 16)

meta <- read.delim('02-meta.txt',head=F,stringsAsFactors=F)
topdim <- 1:30 

seuratObj.list <- lapply(meta$V1, function(x){
  readRDS(sprintf("%s/%s_seurat.rds",indir, x))  
})

# features <- SelectIntegrationFeatures(object.list = seuratObj.list, nfeatures = 3000)
seuratObj.list <- lapply(X = seuratObj.list, FUN = SCTransform, method = "glmGamPoi")
top.features <- SelectIntegrationFeatures(object.list = seuratObj.list, nfeatures = 3000)
seuratObj.list <- PrepSCTIntegration(object.list = seuratObj.list, anchor.features = top.features)
seuratObj.list <- lapply(X = seuratObj.list, FUN = RunPCA, features = top.features)

integrated.anchors <- FindIntegrationAnchors(object.list = seuratObj.list, normalization.method = "SCT", 
    anchor.features = top.features, dims = topdim, reduction = "rpca", k.anchor = 10)
integrated <- IntegrateData(anchorset = integrated.anchors, normalization.method = "SCT", dims = topdim)

#saveRDS(integrated, file="integration_1.rds")
DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
# integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunTSNE(integrated, dims=topdim, dim.embed=3) ##  , 
integrated <- RunUMAP(integrated, dims=topdim, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40) ## n.components=3, dims=topdim, min.dist=0.5, n.neighbors=10, verbose=FALSE, umap.method="umap-learn"
integrated <- FindNeighbors(integrated, dims=topdim, reduction="pca") ## , dims=topdim, k.param=20, reduction='umap', compute.SNN=TRUE, force.recalc=TRUE
integrated <- FindClusters(integrated) ## , resolution=0.2
#saveRDS(integrated, file="integration_2.rds")

integrated@misc$geneName <- seuratObj.list[[1]]@misc$geneName

# Normalize RNA data for visualization purposes
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)

saveRDS(integrated, file="integration_3.rds")

## Identification of all markers for each cluster
top.markers <- FindAllMarkers(integrated, min.pct=0.1, logfc.threshold=0.1, return.thresh=0.1, only.pos=TRUE) ## return.thresh=0.01, test.use="wilcox" 
top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
top.markers$name <- integrated@misc$geneName[top.markers$gene]
## only keep top 100
topmarkers <- top.markers[top.markers$p_val_adj < 0.05,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt=avg_log2FC)
top20 <- topmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=20, wt=avg_log2FC)

integrated@misc$markerGenes <- top.markers
integrated@misc$topMarker <- topmarkers
integrated@misc$top20Marker <- top20

DefaultAssay(integrated) <- "integrated"
saveRDS(integrated, sprintf("%s/integration_4.rds",outdir))
