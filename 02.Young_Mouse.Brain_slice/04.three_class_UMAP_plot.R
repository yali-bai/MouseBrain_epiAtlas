library(Seurat)
library("glmGamPoi")
library(dplyr)
library(future)
library(presto)
library(stringr)
library(getopt)

exc.obj = readRDS("exc_obj.label_transfer_twice.rds")
inte.col <- c("Our"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#89C75F", "Astro-TE NN" = "#0C727C",
                               "L2/3 IT CTX Glut" = "#ed1941", "Sst Gaba" = "#2585a6", "L6 IT CTX Glut" = "#89288F",
                               "OPC NN" = "#F47D2B", "L5 ET CTX Glut" = "#FEE500", "Microglia NN" = "#f26b85",
                               "Peri NN" = "#90D5E4", "CA2-FC-IG Glut" = "#f173ac", "Lamp5 Gaba" = "#00ae9d",
                               "Endo NN" = "#1d953f", "DG Glut" = "#FF6800", "L4/5 IT CTX Glut" = "#009ad6",
                               "Vip Gaba" = "#DEA0FD", "L6 CT CTX Glut" = "#6E4B9E", "Pvalb Gaba" = "#65c294",
                               "CLA-EPd-CTX Car3 Glut" = "#AA0DFE", "L2/3 IT RSP Glut" = '#e74c3c', "CA3 Glut" = "#00538A",
                               "Ependymal NN"="#A6BDD7","CHOR NN"="#B32851","HPF CR Glut"="#F6768E","CA1-ProS Glut"="#b2d235",
                               "L5 NP CTX Glut"="#F4C800","SMC NN"="#C06CAB","VLMC NN"="#007947","L4 RSP-ACA Glut"="#8A9FD1","DG-PIR Ex IMN"='#5AC2F1',
                               "L5 IT CTX Glut" = "#FFB300", "L6b CTX Glut" = '#3498db',"Lamp5 Lhx6 Gaba"='#057771',"RHP-COA Ndnf Gaba"="#e4c6d0",
                               "Sncg Gaba"='#915ce5',"Sst Chodl Gaba"="#ffc773","STR D2 Gaba"="#88c4e8","OB-STR-CTX Inh IMN"="#eb7f54")

Idents(exc.obj) <- exc.obj$subclass
levels(exc.obj)
pdf(paste0("three_class_UMAP_plot/","exc.integration.subclass.UMAP_color_both.pdf"),width = 15,height = 13)
DimPlot(exc.obj, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
exc.obj@meta.data$subclass_color_our_only_label_transfer = exc.obj$subclass
exc.obj@meta.data$subclass_color_our_only_label_transfer[which(is.na(str_match(rownames(exc.obj@meta.data),"Mouses")))] <- "Zhuang"
Idents(exc.obj) <- exc.obj$subclass_color_our_only_label_transfer
levels(exc.obj)
pdf(paste0("three_class_UMAP_plot/","exc.integration.subclass.UMAP_color_our_only.pdf"),width = 15,height = 13)
DimPlot(exc.obj, group.by = "subclass_color_our_only_label_transfer",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(exc.obj)))
dev.off()


inh.obj = readRDS("inh_obj.label_transfer_twice.rds")
Idents(inh.obj) <- inh.obj$subclass
levels(inh.obj)
pdf(paste0("three_class_UMAP_plot/","inh.integration.subclass.UMAP_color_both.pdf"),width = 15,height = 13)
DimPlot(inh.obj, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()
## color our only ##
inh.obj@meta.data$subclass_color_our_only_label_transfer = inh.obj$subclass
inh.obj@meta.data$subclass_color_our_only_label_transfer[which(is.na(str_match(rownames(inh.obj@meta.data),"Mouses")))] <- "Zhuang"
Idents(inh.obj) <- inh.obj$subclass_color_our_only_label_transfer
levels(inh.obj)
pdf(paste0("three_class_UMAP_plot/","inh.integration.subclass.UMAP_color_our_only.pdf"),width = 15,height = 13)
DimPlot(inh.obj, group.by = "subclass_color_our_only_label_transfer",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(inh.obj)))
dev.off()


non.obj = readRDS("non_obj.label_transfer_twice.rds")
Idents(non.obj) <- non.obj$subclass
levels(non.obj)
pdf(paste0("three_class_UMAP_plot/","non.integration.subclass.UMAP_color_both.pdf"),width = 15,height = 13)
DimPlot(non.obj, group.by = "subclass",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE)
dev.off()

## color our only ##
non.obj@meta.data$subclass_color_our_only_label_transfer = non.obj$subclass
non.obj@meta.data$subclass_color_our_only_label_transfer[which(is.na(str_match(rownames(non.obj@meta.data),"Mouses")))] <- "Zhuang"
Idents(non.obj) <- non.obj$subclass_color_our_only_label_transfer
levels(non.obj)
pdf(paste0("three_class_UMAP_plot/","non.integration.subclass.UMAP_color_our_only.pdf"),width = 15,height = 13)
DimPlot(non.obj, group.by = "subclass_color_our_only_label_transfer",cols = inte.col, label=FALSE, pt.size =3,seed=1100,label.box=T,label.size = 6,raster=FALSE,order=rev(levels(non.obj)))
dev.off()


