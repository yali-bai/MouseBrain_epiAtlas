##### import packages
library(Seurat)
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)

##### data prepare ##### 
##### get three class sample id
## read rds of seurat object which integrated our RNA with Zeng 10x
integrated.obj = readRDS("../01.Young_Mouse/01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/Joint_Cabernet.with_celltype.rds")
exc.idx = rownames(integrated.obj@meta.data[which(!is.na(str_match(integrated.obj$three_class,"Exc"))),])
inh.idx = rownames(integrated.obj@meta.data[which(!is.na(str_match(integrated.obj$three_class,"Inh"))),])
non.idx = setdiff(rownames(integrated.obj@meta.data),c(exc.idx,inh.idx))

##### split our MERFISH data
Joint_Cabernet_slice_RNA.df <- readRDS("brain_slice.raw_count.with_QC_filter.rds")
rownames(Joint_Cabernet_slice_RNA.df) <- unlist(lapply(rownames(Joint_Cabernet_slice_RNA.df), function(x) strsplit(x,"\\.")[[1]][1]))

our_exc.df = Joint_Cabernet_slice_RNA.df[,intersect(unlist(lapply(exc.idx, function(x) strsplit(x,"@@_")[[1]][2])),colnames(Joint_Cabernet_slice_RNA.df))]
our_inh.df = Joint_Cabernet_slice_RNA.df[,intersect(unlist(lapply(inh.idx, function(x) strsplit(x,"@@_")[[1]][2])),colnames(Joint_Cabernet_slice_RNA.df))]
our_non.df = Joint_Cabernet_slice_RNA.df[,intersect(unlist(lapply(non.idx, function(x) strsplit(x,"@@_")[[1]][2])),colnames(Joint_Cabernet_slice_RNA.df))]

##### split zhuang MERFISH data
seurat_subset_region <- readRDS("../03.data/02.metainfo/02.Young_Mouse.Brain_slice/subset.z_axis_located_on_7.33.cortex_and_hippo.rds")
zhuang.df <- as.matrix(seurat_subset_region@assays$RNA@counts)
exc.idx = rownames(seurat_subset_region@meta.data[which(!is.na(str_match(seurat_subset_region$subclass_transfer,"Glut"))),])
exc.idx = c(exc.idx,rownames(seurat_subset_region@meta.data[which(!is.na(str_match(seurat_subset_region$subclass_transfer,"DG-PIR Ex IMN"))),]))
inh.idx = rownames(seurat_subset_region@meta.data[which(!is.na(str_match(seurat_subset_region$subclass_transfer,"Gaba"))),])
inh.idx = c(inh.idx,rownames(seurat_subset_region@meta.data[which(!is.na(str_match(seurat_subset_region$subclass_transfer,"OB-STR-CTX Inh IMN"))),]))
non.idx = setdiff(rownames(seurat_subset_region@meta.data),c(exc.idx,inh.idx))
zhuang_exc.df = zhuang.df[,intersect(exc.idx,colnames(zhuang.df))]
zhuang_inh.df = zhuang.df[,intersect(inh.idx,colnames(zhuang.df))]
zhuang_non.df = zhuang.df[,intersect(non.idx,colnames(zhuang.df))]

##### integrated 
common_gene <- intersect(rownames(seurat_subset_region@assays$RNA@counts),rownames(Joint_Cabernet_slice_RNA.df))

integrate = function(our_matrix,zhuang_matrix,common_gene){
    our_matrix = our_matrix[common_gene,]
    zhuang_matrix = zhuang_matrix[common_gene,]

    merged.df <- cbind(zhuang_matrix,our_matrix)
    merged.seuratobj <- CreateSeuratObject(merged.df)
    merged.seuratobj$source <- c(rep("Zhuang",dim(zhuang_matrix)[2]),rep("Joint_Cabernet",dim(our_matrix)[2]))

    ## step 1 ##
    merged.seuratobj.list <- SplitObject(merged.seuratobj, split.by = "source")
    merged.seuratobj.list <- lapply(X = merged.seuratobj.list, FUN = SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(object.list = merged.seuratobj.list)
    merged.seuratobj.list <- PrepSCTIntegration(object.list = merged.seuratobj.list, anchor.features = features)
    merged.seuratobj.list <- lapply(X = merged.seuratobj.list, FUN = RunPCA, features = features)
    merged.seuratobj.anchors <- FindIntegrationAnchors(object.list = merged.seuratobj.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
    merged.seuratobj.sct <- IntegrateData(anchorset = merged.seuratobj.anchors, normalization.method = "SCT", dims = 1:30)

    ## step 2 ##
    DefaultAssay(merged.seuratobj.sct) <- "integrated"
    merged.seuratobj.sct <- RunPCA(merged.seuratobj.sct, verbose = FALSE)
    merged.seuratobj.sct <- RunUMAP(merged.seuratobj.sct, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.5, n.neighbors=40)
    dis=25
    merged.seuratobj.sct <- RunTSNE(merged.seuratobj.sct, dims=1:30, perplexity=dis)
    merged.seuratobj.sct <- FindNeighbors(merged.seuratobj.sct, dims=1:30, reduction="pca")
    merged.seuratobj.sct <- FindClusters(merged.seuratobj.sct)

    ## step 3 ##
    DefaultAssay(merged.seuratobj.sct) <- "RNA"
    merged.seuratobj.sct <- NormalizeData(merged.seuratobj.sct, verbose = FALSE)

    ## step 4 ##
    merged.seuratobj.sct <- JoinLayers(merged.seuratobj.sct)
    top.markers <- FindAllMarkers(merged.seuratobj.sct, min.pct=0.1, logfc.threshold=0.1, return.thresh=0.1, only.pos=TRUE) ## return.thresh=0.01, test.use="wilcox"
    top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
    top.markers$name <- merged.seuratobj.sct@misc$geneName[top.markers$gene]
    ## only keep top 100
    topmarkers <- top.markers[top.markers$p_val_adj < 0.05,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt=avg_log2FC)
    top20 <- topmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=20, wt=avg_log2FC)

    merged.seuratobj.sct@misc$markerGenes <- top.markers
    merged.seuratobj.sct@misc$topMarker <- topmarkers
    merged.seuratobj.sct@misc$top20Marker <- top20

    ##### class
    merged.seuratobj.sct$cell_type <- "NA"
    merged.seuratobj.sct$cell_type[na.omit(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data)))] <- as.character(seurat_subset_region$cell_type)[-which(is.na(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data))))]

    result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
    Idents(merged.seuratobj.sct) <- merged.seuratobj.sct@meta.data$seurat_clusters
    for (ident_value in levels(merged.seuratobj.sct)) {
        cell_indices <- WhichCells(merged.seuratobj.sct, idents = ident_value)
        if(length(which(is.na(str_match(cell_indices,"TSO")))) == 0){
            next
        }
        table_values <- table(merged.seuratobj.sct@meta.data[cell_indices, "cell_type"])
        if(length(which(names(table_values) == "NA")) > 0){
            table_values = table_values[-which(names(table_values) == "NA")]
       }
       max_variable <- names(table_values)[which.max(table_values)]
       result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
    }

    merged.seuratobj.sct@meta.data$major_celltype = NA
    for(i in 1:nrow(result_df)){
        merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$seurat_clusters == result_df$idents[i]),'major_celltype'] <- result_df$max_variable[i]
    }

    ##### subclass 
    merged.seuratobj.sct$subclass_transfer <- NA
    merged.seuratobj.sct$subclass_transfer[na.omit(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data)))] <- as.character(seurat_subset_region$subclass_transfer)[-which(is.na(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data))))]
    result_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
    Idents(merged.seuratobj.sct) <- merged.seuratobj.sct@meta.data$seurat_clusters
    for (ident_value in levels(merged.seuratobj.sct)) {
        cell_indices <- WhichCells(merged.seuratobj.sct, idents = ident_value)
        if(length(which(is.na(str_match(cell_indices,"TSO")))) == 0){
            next
        }
        table_values <- table(merged.seuratobj.sct@meta.data[cell_indices, "subclass_transfer"])
        if(length(which(names(table_values) == "NA")) > 0){
            table_values = table_values[-which(names(table_values) == "NA")]
        }
        max_variable <- names(table_values)[which.max(table_values)]
        result_df <- rbind(result_df, data.frame(idents = ident_value, max_variable = max_variable))
    }

    ## all ##
    merged.seuratobj.sct@meta.data$subclass = NA
    for(i in 1:nrow(result_df)){
        merged.seuratobj.sct@meta.data[which(merged.seuratobj.sct@meta.data$seurat_clusters == result_df$idents[i]),'subclass'] <- result_df$max_variable[i]
    }
    merged.seuratobj.sct@meta.data$major_brain_region = unlist(lapply(rownames(merged.seuratobj.sct@meta.data), function(x) strsplit(x,"_")[[1]][3]))
    merged.seuratobj.sct@meta.data$major_brain_region[na.omit(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data)))] <- as.character(seurat_subset_region@meta.data$major_brain_region)[-which(is.na(match(rownames(seurat_subset_region@meta.data),rownames(merged.seuratobj.sct@meta.data))))]
    merged.seuratobj.sct@meta.data$major_brain_region_v2 = merged.seuratobj.sct@meta.data$major_brain_region
    merged.seuratobj.sct@meta.data$major_brain_region_v2[which(merged.seuratobj.sct@meta.data$major_brain_region_v2 == "LeftCortex")] = "Isocortex"
    merged.seuratobj.sct@meta.data$major_brain_region_v2[which(merged.seuratobj.sct@meta.data$major_brain_region_v2 == "RightCortex")] = "Isocortex"
    merged.seuratobj.sct@meta.data$major_brain_region_v2[which(merged.seuratobj.sct@meta.data$major_brain_region_v2 == "LeftHippo")] = "Hippocampus"
    merged.seuratobj.sct@meta.data$major_brain_region_v2[which(merged.seuratobj.sct@meta.data$major_brain_region_v2 == "RightHippo")] = "Hippocampus"
    
    return(merged.seuratobj.sct)
}


exc.obj = integrate(our_exc.df,zhuang_exc.df,common_gene)
inh.obj = integrate(our_inh.df,zhuang_inh.df,common_gene)
non.obj = integrate(our_non.df,zhuang_non.df,common_gene)

saveRDS(exc.obj,file="exc_obj.label_transfer_twice.rds")
saveRDS(inh.obj,file="inh_obj.label_transfer_twice.rds")
saveRDS(non.obj,file="non_obj.label_transfer_twice.rds")