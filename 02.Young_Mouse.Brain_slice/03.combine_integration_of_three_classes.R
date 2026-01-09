
##### 01. import packages #####
library(Seurat)
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)
library(ggplot2)

##### integrated RNA matrix
exc.obj = readRDS("exc_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
inh.obj = readRDS("inh_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
non.obj = readRDS("non_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
exc.matrix = as.matrix(exc.obj@assays$RNA$counts)
inh.matrix = as.matrix(inh.obj@assays$RNA$counts)
non.matrix = as.matrix(non.obj@assays$RNA$counts)

merged.df <- cbind(exc.matrix,inh.matrix,non.matrix)
merged.seuratobj <- CreateSeuratObject(merged.df)
merged.seuratobj$source = "Zhuang"
merged.seuratobj$source[which(!is.na(str_match(rownames(merged.seuratobj@meta.data),"TSO")))] = "Joint_Cabernet"

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
top.markers <- FindAllMarkers(merged.seuratobj.sct, min.pct=0.1, logfc.threshold=0.1, return.thresh=0.1, only.pos=TRUE) 
top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
top.markers$name <- merged.seuratobj.sct@misc$geneName[top.markers$gene]
## only keep top 100
topmarkers <- top.markers[top.markers$p_val_adj < 0.05,] %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt=avg_log2FC)
top20 <- topmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=20, wt=avg_log2FC)

merged.seuratobj.sct@misc$markerGenes <- top.markers
merged.seuratobj.sct@misc$topMarker <- topmarkers
merged.seuratobj.sct@misc$top20Marker <- top20

saveRDS(merged.seuratobj.sct,file="integrated.rds")


##### add sample info 
exc.obj = readRDS("exc_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
inh.obj = readRDS("inh_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
non.obj = readRDS("non_obj.label_transfer_twice.loci_transfer.the_nearst_1_cell.rds")
result.df = exc.obj@meta.data
result.df = rbind(result.df,inh.obj@meta.data)
result.df = rbind(result.df,non.obj@meta.data)

result.df = result.df[rownames(merged.seuratobj.sct@meta.data),]
merged.seuratobj.sct@meta.data = cbind(merged.seuratobj.sct@meta.data,result.df[,9:dim(result.df)[2]])

if (!dir.exists(paste0("cluster_plot"))) {  
    # 如果文件夹不存在，则创建它  
    dir.create(paste0("cluster_plot"))  
}

### 06. plot ###
## 01.  integrated seurat cluster ##
pdf("./cluster_plot/integrated.seurat_clusters.tsne.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "seurat_clusters",label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
pdf("./cluster_plot/integrated.seurat_clusters.UMAP.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "seurat_clusters",label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## 02. major class ##
## all ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype
levels(merged.seuratobj.sct)

inte.col <- c("Zhuang"='lightgrey',"Joint_Cabernet"='lightgrey',"endothelial cell"="#f173ac","oligodendrocyte"="#8A9FD1","GABAergic neuron"="#ed1941","microglial cell"="#2585a6","astrocyte"="#89288F","oligodendrocyte precursor cell"="#FF6800","glutamatergic neuron"="#90D5E4","pericyte"="#da765b","ependymal cell"="#00ae9d","choroid plexus epithelial cell"="#3283FE","smooth muscle cell"="#8552a1","vascular leptomeningeal cell"="#00538A","neuroblast (sensu Vertebrata)"="#FEE500")
pdf("./cluster_plot/integration.major_celltype.tsne_all.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_celltype",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
pdf("./cluster_plot/integration.major_celltype.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
merged.seuratobj.sct@meta.data$major_celltype_color_our_only = merged.seuratobj.sct$major_celltype
merged.seuratobj.sct@meta.data$major_celltype_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype_color_our_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.major_celltype.tsne_color_our_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_celltype_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf("./cluster_plot/integration.major_celltype.UMAP_color_our_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
merged.seuratobj.sct@meta.data$major_celltype_color_zhuang_only = merged.seuratobj.sct$major_celltype
merged.seuratobj.sct@meta.data$major_celltype_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype_color_zhuang_only
levels(merged.seuratobj.sct)

pdf("./cluster_plot/integration.major_celltype.tsne_color_zhuang_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_celltype_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()
pdf("./cluster_plot/integration.major_celltype.UMAP_color_zhuang_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

## brain region ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2
levels(merged.seuratobj.sct)

inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Isocortex"="#3497dc", "Hippocampus"="#2bc96d", "LeftCortex"="#8A9FD1", "LeftHippo"="#ed1941", "RightCortex"="#FEE500", "RightHippo"="#b2d235")

pdf("./cluster_plot/integration.major_brain_region_v2.tsne_all.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf("./cluster_plot/integration.major_brain_region_v2.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color our only ##
merged.seuratobj.sct@meta.data$major_brain_region_v2_color_our_only <- merged.seuratobj.sct@meta.data$major_brain_region_v2
merged.seuratobj.sct@meta.data$major_brain_region_v2_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2_color_our_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.major_brain_region_v2.tsne_color_our_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf("./cluster_plot/integration.major_brain_region_v2.UMAP_color_our_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
merged.seuratobj.sct@meta.data$major_brain_region_v2_color_zhuang_only <- merged.seuratobj.sct@meta.data$major_brain_region_v2
merged.seuratobj.sct@meta.data$major_brain_region_v2_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2_color_zhuang_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.major_brain_region_v2.tsne_color_zhuang_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Hippocampus","Isocortex","Joint_Cabernet"))
dev.off()
pdf("./cluster_plot/integration.major_brain_region_v2.UMAP_color_zhuang_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Hippocampus","Isocortex","Joint_Cabernet"))
dev.off()

## subclass ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass
levels(merged.seuratobj.sct)

inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#f173ac", "Astro-TE NN" = "#0C727C",
                 "L2/3 IT CTX Glut" = "#ed1941", "Sst Gaba" = "#2585a6", "L6 IT CTX Glut" = "#89288F", 
                 "OPC NN" = "#F47D2B", "L5 ET CTX Glut" = "#FEE500", "Microglia NN" = "#da765b",  
                 "Peri NN" = "#90D5E4", "CA2-FC-IG Glut" = "#00538A", "Lamp5 Gaba" = "#00ae9d", 
                 "Endo NN" = "#1d953f", "DG Glut" = "#FF6800", "L4/5 IT CTX Glut" = "#009ad6", 
                 "Vip Gaba" = "#DEA0FD", "L6 CT CTX Glut" = "#6E4B9E", "Pvalb Gaba" = "#65c294", 
                 "CLA-EPd-CTX Car3 Glut" = "#AA0DFE", "L2/3 IT RSP Glut" = '#e74c3c', "CA3 Glut" = "#89C75F",
                 "Ependymal NN"="#A6BDD7","CHOR NN"="#B32851","HPF CR Glut"="#F6768E","CA1-ProS Glut"="#b2d235",
                 "L5 NP CTX Glut"="#F4C800","SMC NN"="#C06CAB","VLMC NN"="#007947","L4 RSP-ACA Glut"="#8A9FD1","DG-PIR Ex IMN"='#5AC2F1FF',
                 "L5 IT CTX Glut" = "#FFB300", "L6b CTX Glut" = '#3498db',"Lamp5 Lhx6 Gaba"='#00a6ac',"RHP-COA Ndnf Gaba"="#e4c6d0",
                 "Sncg Gaba"="#f9906f","Sst Chodl Gaba"="#ffc773","STR D2 Gaba"="#88c4e8","OB-STR-CTX Inh IMN"="#eb7f54")


pdf("./cluster_plot/integration.subclass.tsne_all.pdf",width = 17,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
pdf("./cluster_plot/integration.subclass.UMAP_all.pdf",width = 17,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
merged.seuratobj.sct$subclass_color_our_only <- merged.seuratobj.sct@meta.data$subclass
merged.seuratobj.sct$subclass_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"

Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass_color_our_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.subclass.tsne_color_our_only.pdf",width = 17,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()
pdf("./cluster_plot/integration.subclass.UMAP_color_our_only.pdf",width = 17,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
merged.seuratobj.sct$subclass_color_zhuang_only <- merged.seuratobj.sct@meta.data$subclass
merged.seuratobj.sct$subclass_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"

Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass_color_zhuang_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.subclass.tsne_color_zhuang_only.pdf",width = 17,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "subclass_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()
pdf("./cluster_plot/integration.subclass.UMAP_color_zhuang_only.pdf",width = 17,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

## source ##
merged.seuratobj.sct$source = NA
merged.seuratobj.sct$source[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] = "Zhuang"
merged.seuratobj.sct$source[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] = "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$source
levels(merged.seuratobj.sct)

inte.col <- setNames(c('#3498db','#e74c3c'),
                    c("Zhuang","Joint_Cabernet"))
pdf("./cluster_plot/integration.source.tsne_all.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "source",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang"))
dev.off()
pdf("./cluster_plot/integration.source.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "source",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang"))
dev.off()

## 07. neuron vs non-neuron ##
merged.seuratobj.sct$neuron_type = "Non_neuron"
merged.seuratobj.sct$neuron_type[which(!is.na(str_match(merged.seuratobj.sct$subclass,"Gaba")))] <- "Neuron"
merged.seuratobj.sct$neuron_type[which(!is.na(str_match(merged.seuratobj.sct$subclass,"Glut")))] <- "Neuron"
merged.seuratobj.sct$neuron_type[which(!is.na(str_match(merged.seuratobj.sct$subclass,"DG-PIR Ex IMN")))] <- "Neuron"
merged.seuratobj.sct$neuron_type[which(!is.na(str_match(merged.seuratobj.sct$subclass,"OB-STR-CTX Inh IMN")))] <- "Neuron"
## plot ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type
levels(merged.seuratobj.sct)


inte.col <- c("Neuron" = "#19bc9b","Non_neuron" = "#f39b11","Zhuang"= "lightgrey","Joint_Cabernet"= "lightgrey")

pdf("./cluster_plot/integration.neuron_vs_nonneuron.tsne_all.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "neuron_type",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
pdf("./cluster_plot/integration.neuron_vs_nonneuron.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
merged.seuratobj.sct$neuron_type_color_our_only <- merged.seuratobj.sct$neuron_type
merged.seuratobj.sct$neuron_type_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type_color_our_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.neuron_vs_nonneuron.tsne_color_our_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "neuron_type_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Zhuang"))
dev.off()
pdf("./cluster_plot/integration.neuron_vs_nonneuron.UMAP_color_our_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Zhuang"))
dev.off()

## color zhuang only ##
merged.seuratobj.sct$neuron_type_color_zhuang_only <- merged.seuratobj.sct$neuron_type
merged.seuratobj.sct$neuron_type_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type_color_zhuang_only
levels(merged.seuratobj.sct)


pdf("./cluster_plot/integration.neuron_vs_nonneuron.tsne_color_zhuang_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "neuron_type_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Joint_Cabernet"))
dev.off()
pdf("./cluster_plot/integration.neuron_vs_nonneuron.UMAP_color_zhuang_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Joint_Cabernet"))
dev.off()

##### 08.three class #####
inte.col <- c("Exc"='#e74c3c',"Inh" = '#33a3dc',"Non_neuron"="orange","Zhuang"='lightgray',"Joint_Cabernet"='lightgray')

merged.seuratobj.sct$three_class = "Non_neuron"
merged.seuratobj.sct$three_class[which(!is.na(str_match(merged.seuratobj.sct$subclass,"Glut")))] = "Exc"
merged.seuratobj.sct$three_class[which(!is.na(str_match(merged.seuratobj.sct$subclass,"Gaba")))] = "Inh"
merged.seuratobj.sct$three_class[which(merged.seuratobj.sct$subclass=="DG-PIR Ex IMN")] <- "Exc"
merged.seuratobj.sct$three_class[which(merged.seuratobj.sct$subclass=="OB-STR-CTX Inh IMN")] <- "Inh"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class
levels(merged.seuratobj.sct)

pdf("./cluster_plot/integration.three_class.tsne_all.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "three_class",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
pdf("./cluster_plot/integration.three_class.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

##### color our only
merged.seuratobj.sct$three_class_color_our_only <- merged.seuratobj.sct$three_class
merged.seuratobj.sct$three_class_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class_color_our_only
levels(merged.seuratobj.sct)

pdf("./cluster_plot/integration.three_class.tsne_color_our_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "three_class_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = rev(levels(merged.seuratobj.sct)))
dev.off()
pdf("./cluster_plot/integration.three_class.UMAP_color_our_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class_color_our_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
merged.seuratobj.sct$three_class_color_zhuang_only <- merged.seuratobj.sct$three_class
merged.seuratobj.sct$three_class_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class_color_zhuang_only
levels(merged.seuratobj.sct)

pdf("./cluster_plot/integration.three_class.tsne_color_zhuang_only.pdf",width = 15,height = 13)
TSNEPlot(merged.seuratobj.sct, group.by = "three_class_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = levels(merged.seuratobj.sct))
dev.off()
pdf("./cluster_plot/integration.three_class.UMAP_color_zhuang_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = levels(merged.seuratobj.sct))
dev.off()

saveRDS(merged.seuratobj.sct,file="integrated.with_sampleinfo.rds")

##### map
exc.loci = data.frame(Embeddings(exc.obj, reduction = 'umap')[Cells(exc.obj),])
inh.loci = data.frame(Embeddings(inh.obj, reduction = 'umap')[Cells(inh.obj),])
non.loci = data.frame(Embeddings(non.obj, reduction = 'umap')[Cells(non.obj),])
loci_inte.df = rbind(exc.loci,inh.loci,non.loci)
loci_inte.df = loci_inte.df[rownames(merged.seuratobj.sct@meta.data),]
loci_inte.df = as.matrix(loci_inte.df)
merged.seuratobj.sct@reductions$umap@cell.embeddings <- loci_inte.df
saveRDS(merged.seuratobj.sct,file="merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")

##### plot #####
if (!dir.exists(paste0("spatial_plot_",1))) {  
    # 如果文件夹不存在，则创建它  
    dir.create(paste0("spatial_plot_",1))  
}

## 01.  integrated seurat cluster ##
pdf(paste0("./spatial_plot_",1,"/integrated.seurat_clusters.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "seurat_clusters",label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## 02. neuron vs non-neuron ##
## plot ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type
levels(merged.seuratobj.sct)

our.col <- setNames(c('#3498db','#e74c3c'),
                    c("Neuron","Non-neuron"))

pdf(paste0("./spatial_plot_",1,"/integration.neuron_vs_nonneuron.umap_all.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type",cols = our.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type_color_our_only
levels(merged.seuratobj.sct)
our.col <- setNames(c('#3498db','#e74c3c','lightgray'),
                    c("Neuron","Non-neuron","Zhuang"))

pdf(paste0("./spatial_plot_",1,"/integration.neuron_vs_nonneuron.umap_color_our_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type_color_our_only",cols = our.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Zhuang"))
dev.off()

## color zhuang only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$neuron_type_color_zhuang_only
levels(merged.seuratobj.sct)
our.col <- setNames(c('#3498db','#e74c3c','lightgray'),
                    c("Neuron","Non-neuron","Joint_Cabernet"))

pdf(paste0("./spatial_plot_",1,"/integration.neuron_vs_nonneuron.umap_color_zhuang_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "neuron_type_color_zhuang_only",cols = our.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = c("Neuron","Non-neuron","Joint_Cabernet"))
dev.off()

## 03. cell type ##
inte.col <- c("Zhuang"='lightgrey',"Joint_Cabernet"='lightgrey',"endothelial cell"="#f173ac","oligodendrocyte"="#8A9FD1","GABAergic neuron"="#ed1941","microglial cell"="#2585a6","astrocyte"="#89288F","oligodendrocyte precursor cell"="#FF6800","glutamatergic neuron"="#90D5E4","pericyte"="#da765b","ependymal cell"="#00ae9d","choroid plexus epithelial cell"="#3283FE","smooth muscle cell"="#8552a1","vascular leptomeningeal cell"="#00538A","neuroblast (sensu Vertebrata)"="#FEE500")

## all ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_celltype.umap_all.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype_color_our_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_celltype.umap_color_our_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype_color_our_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype_color_zhuang_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_celltype.umap_color_zhuang_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_celltype_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

##### 05. major brain region #####
inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Isocortex"="#C06CAB", "Hippocampus"="#90D5E4", "LeftCortex"="#8A9FD1", "LeftHippo"="#ed1941", "RightCortex"="#FEE500", "RightHippo"="#b2d235")

## all: Isocortex Hippocampus ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_brain_region_v2.umap_all.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color our only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2_color_our_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_brain_region_v2.umap_color_our_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_our_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_brain_region_v2_color_zhuang_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.major_brain_region_v2.umap_color_zhuang_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "major_brain_region_v2_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Fiber_tracts","Hypothalamus","Olfactory","Striatum","Ventricular_systems","Cortical_subplate","Pallidum","Thalamus","Isocortex","Hippocampus","Joint_Cabernet","NA"))
dev.off()

##### 06. subclass_transfer #####
inte.col <- c("Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#f173ac", "Astro-TE NN" = "#0C727C",
                 "L2/3 IT CTX Glut" = "#ed1941", "Sst Gaba" = "#2585a6", "L6 IT CTX Glut" = "#89288F", 
                 "OPC NN" = "#F47D2B", "L5 ET CTX Glut" = "#FEE500", "Microglia NN" = "#da765b",  
                 "Peri NN" = "#90D5E4", "CA2-FC-IG Glut" = "#00538A", "Lamp5 Gaba" = "#00ae9d", 
                 "Endo NN" = "#1d953f", "DG Glut" = "#FF6800", "L4/5 IT CTX Glut" = "#009ad6", 
                 "Vip Gaba" = "#DEA0FD", "L6 CT CTX Glut" = "#6E4B9E", "Pvalb Gaba" = "#65c294", 
                 "CLA-EPd-CTX Car3 Glut" = "#AA0DFE", "L2/3 IT RSP Glut" = '#e74c3c', "CA3 Glut" = "#89C75F",
                 "Ependymal NN"="#A6BDD7","CHOR NN"="#B32851","HPF CR Glut"="#F6768E","CA1-ProS Glut"="#b2d235",
                 "L5 NP CTX Glut"="#F4C800","SMC NN"="#C06CAB","VLMC NN"="#007947","L4 RSP-ACA Glut"="#8A9FD1","DG-PIR Ex IMN"='#5AC2F1FF',
                 "L5 IT CTX Glut" = "#FFB300", "L6b CTX Glut" = '#3498db',"Lamp5 Lhx6 Gaba"='#00a6ac',"RHP-COA Ndnf Gaba"="#e4c6d0",
                 "Sncg Gaba"="#f9906f","Sst Chodl Gaba"="#ffc773","STR D2 Gaba"="#88c4e8","OB-STR-CTX Inh IMN"="#eb7f54")

## all ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.subclass.umap_all.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass_color_our_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.subclass.umap_color_our_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_color_our_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$subclass_color_zhuang_only
levels(merged.seuratobj.sct)

pdf(paste0("./spatial_plot_",1,"/integration.subclass.umap_color_zhuang_only.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "subclass_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=levels(merged.seuratobj.sct))
dev.off()

##### 07. source #####
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$source
levels(merged.seuratobj.sct)

inte.col <- setNames(c('#3498db','#e74c3c'),
                    c("Zhuang","Joint_Cabernet"))
pdf(paste0("./spatial_plot_",1,"/integration.source.umap_all.pdf"),width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "source",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang"))
dev.off()

##### 08.three class #####
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class
levels(merged.seuratobj.sct)

inte.col <- setNames(c("#FFB300",'#3498db','#e74c3c','lightgray','lightgray'),
                    c("Inh","Exc","Non-neuron","Zhuang","Joint_Cabernet"))

pdf("./spatial_plot_1/integration.three_class.UMAP_all.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

##### color our only
#merged.seuratobj.sct$three_class_color_our_only <- merged.seuratobj.sct$three_class
#merged.seuratobj.sct$three_class_color_our_only[which(is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Zhuang"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class_color_our_only
levels(merged.seuratobj.sct)

pdf("./spatial_plot_1/integration.three_class.UMAP_color_our_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class_color_our_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = rev(levels(merged.seuratobj.sct)))
dev.off()

## color zhuang only ##
#merged.seuratobj.sct$three_class_color_zhuang_only <- merged.seuratobj.sct$three_class
#merged.seuratobj.sct$three_class_color_zhuang_only[which(!is.na(str_match(rownames(merged.seuratobj.sct@meta.data),"TSO")))] <- "Joint_Cabernet"
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$three_class_color_zhuang_only
levels(merged.seuratobj.sct)

pdf("./spatial_plot_1/integration.three_class.UMAP_color_zhuang_only.pdf",width = 15,height = 13)
DimPlot(merged.seuratobj.sct, group.by = "three_class_color_zhuang_only",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order = levels(merged.seuratobj.sct))
dev.off()


##### color only one cluster #####
Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$major_celltype
levels(merged.seuratobj.sct)

color_one_cluster_plot <- function(cl, data){
    if (!dir.exists(paste0("color_one_cluster_plot_",1))) {  
        # 如果文件夹不存在，则创建它  
        dir.create(paste0("color_one_cluster_plot_",1))  
    }
    if(data == "all"){
        merged.seuratobj.sct$color_one_cluster <- "other_cluster"
        merged.seuratobj.sct$color_one_cluster[intersect(which(merged.seuratobj.sct$major_celltype == cl),which(merged.seuratobj.sct$source == "Zhuang"))] <- "Zhuang"
        merged.seuratobj.sct$color_one_cluster[intersect(which(merged.seuratobj.sct$major_celltype == cl),which(merged.seuratobj.sct$source == "Joint_Cabernet"))] <- "Joint_Cabernet"
        Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$color_one_cluster
        levels(merged.seuratobj.sct)
        inte.col <- c("other_cluster"='lightgrey',"Zhuang"="#415284", "Joint_Cabernet"="#EE934E") 
        pic <- DimPlot(merged.seuratobj.sct, group.by = "color_one_cluster",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang","other_cluster"))
        ggsave(pic,file=paste0("./color_one_cluster_plot_",1,"/integration.only_color_",cl,".umap_",data,".color_by_source.pdf"),height=13,width=15)

        merged.seuratobj.sct$color_one_cluster <- "other_cluster"
        merged.seuratobj.sct$color_one_cluster[which(merged.seuratobj.sct$major_celltype == cl)] <- cl
    }
    else{
        merged.seuratobj.sct$color_one_cluster <- "other_cluster"
        merged.seuratobj.sct$color_one_cluster[intersect(which(merged.seuratobj.sct$major_celltype == cl),which(merged.seuratobj.sct$source == data))] <- cl
    }
    Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$color_one_cluster
    levels(merged.seuratobj.sct)

    inte.col <- c("other_cluster"='lightgrey',"endothelial cell"="#f173ac","oligodendrocyte"="#8A9FD1","GABAergic neuron"="#ed1941","microglial cell"="#2585a6","astrocyte"="#89288F","oligodendrocyte precursor cell"="#FF6800","glutamatergic neuron"="#90D5E4","pericyte"="#da765b","ependymal cell"="#00ae9d","choroid plexus epithelial cell"="#3283FE","smooth muscle cell"="#8552a1","vascular leptomeningeal cell"="#00538A","neuroblast (sensu Vertebrata)"="#FEE500")

    pic <- DimPlot(merged.seuratobj.sct, group.by = "color_one_cluster",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c(cl,"other_cluster"))
    ggsave(pic,file=paste0("./color_one_cluster_plot_",1,"/integration.only_color_",cl,".umap_",data,".pdf"),height=13,width=15)
}

for (cl in unique(merged.seuratobj.sct$major_celltype)){
    for(data in c("all","Zhuang","Joint_Cabernet")){
        color_one_cluster_plot(cl,data)
    }
}

##### color only one subclass #####
color_one_subclass_plot <- function(cl, data){
    if (!dir.exists(paste0("color_one_subclass_plot_",1))) {  
        # creat directory if it does not exist  
        dir.create(paste0("color_one_subclass_plot_",1))  
    }
    if(data == "all"){
        merged.seuratobj.sct$color_one_subclass <- "other_subclass"
        merged.seuratobj.sct$color_one_subclass[intersect(which(merged.seuratobj.sct$subclass == cl),which(merged.seuratobj.sct$source == "Zhuang"))] <- "Zhuang"
        merged.seuratobj.sct$color_one_subclass[intersect(which(merged.seuratobj.sct$subclass == cl),which(merged.seuratobj.sct$source == "Joint_Cabernet"))] <- "Joint_Cabernet"
        Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$color_one_subclass
        levels(merged.seuratobj.sct)
        inte.col <- c("other_subclass"='lightgrey',"Zhuang"="#415284", "Joint_Cabernet"="#EE934E") 
        pic <- DimPlot(merged.seuratobj.sct, group.by = "color_one_subclass",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang","other_cluster"))
        ggsave(pic,file=paste0("./color_one_subclass_plot_",1,"/integration.only_color_",str_replace(cl, "/", " "),".umap_",data,".color_by_source.pdf"),height=13,width=15)

        merged.seuratobj.sct$color_one_subclass <- "other_subclass"
        merged.seuratobj.sct$color_one_subclass[which(merged.seuratobj.sct$subclass == cl)] <- cl
    }
    else{
        merged.seuratobj.sct$color_one_subclass <- "other_subclass"
        merged.seuratobj.sct$color_one_subclass[intersect(which(merged.seuratobj.sct$subclass == cl),which(merged.seuratobj.sct$source == data))] <- cl
    }
    Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$color_one_subclass
    levels(merged.seuratobj.sct)

    inte.col <- c("other_subclass"='lightgrey',"Joint_Cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#f173ac", "Astro-TE NN" = "#0C727C",
                 "L2/3 IT CTX Glut" = "#ed1941", "Sst Gaba" = "#2585a6", "L6 IT CTX Glut" = "#89288F", 
                 "OPC NN" = "#F47D2B", "L5 ET CTX Glut" = "#FEE500", "Microglia NN" = "#da765b",  
                 "Peri NN" = "#90D5E4", "CA2-FC-IG Glut" = "#00538A", "Lamp5 Gaba" = "#00ae9d", 
                 "Endo NN" = "#1d953f", "DG Glut" = "#FF6800", "L4/5 IT CTX Glut" = "#009ad6", 
                 "Vip Gaba" = "#DEA0FD", "L6 CT CTX Glut" = "#6E4B9E", "Pvalb Gaba" = "#65c294", 
                 "CLA-EPd-CTX Car3 Glut" = "#AA0DFE", "L2/3 IT RSP Glut" = '#e74c3c', "CA3 Glut" = "#89C75F",
                 "Ependymal NN"="#A6BDD7","CHOR NN"="#B32851","HPF CR Glut"="#F6768E","CA1-ProS Glut"="#b2d235",
                 "L5 NP CTX Glut"="#F4C800","SMC NN"="#C06CAB","VLMC NN"="#007947","L4 RSP-ACA Glut"="#8A9FD1","DG-PIR Ex IMN"='#5AC2F1FF',
                 "L5 IT CTX Glut" = "#FFB300", "L6b CTX Glut" = '#3498db',"Lamp5 Lhx6 Gaba"='#00a6ac',"RHP-COA Ndnf Gaba"="#e4c6d0",
                 "Sncg Gaba"="#f9906f","Sst Chodl Gaba"="#ffc773","STR D2 Gaba"="#88c4e8","OB-STR-CTX Inh IMN"="#eb7f54")

    pic <- DimPlot(merged.seuratobj.sct, group.by = "color_one_subclass",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c(cl,"other_cluster"))
    ggsave(pic,file=paste0("./color_one_subclass_plot_",1,"/integration.only_color_",str_replace(cl, "/", " "),".umap_",data,".pdf"),height=13,width=15)

}

for (cl in unique(merged.seuratobj.sct$subclass)){
    for(data in c("all","Zhuang","Joint_Cabernet")){
        color_one_subclass_plot(cl,data)
    }
}

##### color by source #####
color_by_source_plot <- function(type){  ## "neuron_type"
    if (!dir.exists(paste0("color_by_source_plot_",1))) {  
        # creat directory if it does not exist 
        dir.create(paste0("color_by_source_plot_",1))  
    }
    for (cl in unique(merged.seuratobj.sct@meta.data[,type])){
        merged.seuratobj.sct$color_by_source <- "other_type"
        merged.seuratobj.sct$color_by_source[intersect(which(merged.seuratobj.sct@meta.data[,type] == cl),which(merged.seuratobj.sct$source == "Zhuang"))] <- "Zhuang"
        merged.seuratobj.sct$color_by_source[intersect(which(merged.seuratobj.sct@meta.data[,type] == cl),which(merged.seuratobj.sct$source == "Joint_Cabernet"))] <- "Joint_Cabernet"
        Idents(merged.seuratobj.sct) <- merged.seuratobj.sct$color_by_source
        levels(merged.seuratobj.sct)
        inte.col <- c("other_type"='lightgrey', "Zhuang"="#415284", "Joint_Cabernet"="#EE934E")#"Zhuang"='#3498db',"Joint_Cabernet"='#e74c3c')
        pic <- DimPlot(merged.seuratobj.sct, group.by = "color_by_source",cols = inte.col, label=FALSE, pt.size =1.5,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=c("Joint_Cabernet","Zhuang","other_type"))
        ggsave(pic,file=paste0("./color_by_source_plot_",1,"/integration.only_color_",str_replace(cl, "/", " "),".umap.color_by_source.pdf"),height=13,width=15)
    }
}

for (type in c("neuron_type","major_brain_region_v2")){
    color_by_source_plot(type)
}



