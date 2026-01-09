library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
library(data.table)
library(paletteer)
library(ComplexHeatmap)

seuratObj <- readRDS("integration_4.rds")

##### seurat cluster #####
## plot seurat cluster ##
pdf("seurat_clusters.UMAP.without_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "seurat_clusters",label=FALSE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

pdf("seurat_clusters.UMAP.with_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "seurat_clusters",label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

seuratObj@meta.data$region <- rownames(seuratObj@meta.data)
seuratObj@meta.data$region <- gsub("@@.*","",seuratObj@meta.data$region)

##### data source #####
seuratObj@meta.data$group = seuratObj@meta.data$region
seuratObj@meta.data$group <- ifelse(seuratObj@meta.data$group == "Joint_Cabernet","Joint_Cabernet","Zeng")

## plot data source ##
inte.col=c("Joint_Cabernet" = '#e74c3c',"Zeng"='#33a3dc')
pdf("data_source.UMAP.without_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "group",label=FALSE, cols = inte.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=c("Joint_Cabernet","Zeng"),raster=FALSE)
dev.off()

pdf("data_source.UMAP.with_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "group",label=TRUE, cols = inte.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=c("Joint_Cabernet","Zeng"),raster=FALSE)
dev.off()

# change size
pltd0 <- data.frame(Embeddings(seuratObj, reduction = 'umap')[Cells(seuratObj),],
                    group=ifelse(seuratObj$group == "Joint_Cabernet","Joint_Cabernet","Zeng")) %>% 
    mutate(size=ifelse(group == "Joint_Cabernet", 0.2,0.1)) %>% 
    arrange(group)

pltd0_JC = pltd0[pltd0$group == "Joint_Cabernet",]
pltd0_Z = pltd0[pltd0$group == "Zeng",]
inte.col=c("Joint_Cabernet" = '#e74c3c',"Zeng"="lightgrey")
pdf("data_source.color_Joint_Cabernet_only.UMAP.without_label.change_point_size.pdf",width = 15,height = 13)
ggplot() + 
  geom_point(data=pltd0_Z,aes(x=umap_1,y=umap_2, color=group), size=1) +
  geom_point(data=pltd0_JC,aes(x=umap_1,y=umap_2, color=group), size=3) +  
  scale_color_manual(values = inte.col) + 
  ggpubr::theme_pubr() #+ NoLegend() + NoAxes()
dev.off()

##### region #####
seuratObj@meta.data$region[which(seuratObj@meta.data$group == "Joint_Cabernet")] = unlist(lapply(rownames(seuratObj@meta.data)[which(seuratObj@meta.data$group == "Joint_Cabernet")], function(x) strsplit(x,"_")[[1]][5]))
seuratObj$big_region <- ifelse(seuratObj$region %in% c("Hippo","RightHippo","LeftHippo","HPF"),"Hippo",
                              ifelse(seuratObj$region %in% c("CTX","SuperficialCortex","DeepCortex","LeftCortex","RightCortex"),"Cortex",
                                    ifelse(seuratObj$region %in% c("Olfactory","OLF"),"Olfactory",
                                          ifelse(seuratObj$region %in% c("Striatum","STR"),"Striatum","Cortex"))))

## plot region ##
Idents(seuratObj) <- seuratObj$big_region
levels(seuratObj)
big_region.col = c("Hippo" = "#2bc96d","Cortex"="#3497dc","Olfactory"="#f2c30f","Striatum" = "#e74c3c" ,"Zeng"= "lightgrey","Joint_Cabernet"= "lightgrey")

pdf("big_region.UMAP.without_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region",label=FALSE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=c("Joint_Cabernet","Zeng"),raster=FALSE)
dev.off()

pdf("big_region.UMAP.with_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region",label=TRUE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=c("Joint_Cabernet","Zeng"),raster=FALSE)
dev.off()

seuratObj$big_region_color_Joint_Cabernet_only = seuratObj$big_region
seuratObj$big_region_color_Joint_Cabernet_only[which(seuratObj$group == "Zeng")] = "Zeng"
Idents(seuratObj) <- seuratObj$big_region_color_Joint_Cabernet_only
levels(seuratObj)

pdf("big_region_color_Joint_Cabernet_only.UMAP.without_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region_color_Joint_Cabernet_only",label=FALSE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=rev(levels(seuratObj)),raster=FALSE)
dev.off()

pdf("big_region_color_Joint_Cabernet_only.UMAP.with_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region_color_Joint_Cabernet_only",label=TRUE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=rev(levels(seuratObj)),raster=FALSE)
dev.off()

seuratObj$big_region_color_Zeng_only = seuratObj$big_region
seuratObj$big_region_color_Zeng_only[which(seuratObj$group == "Joint_Cabernet")] = "Joint_Cabernet"
Idents(seuratObj) <- seuratObj$big_region_color_Zeng_only
levels(seuratObj)

pdf("big_region_color_Zeng_only.UMAP.without_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region_color_Zeng_only",label=FALSE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=levels(seuratObj),raster=FALSE)
dev.off()

pdf("big_region_color_Zeng_only.UMAP.with_label.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "big_region_color_Zeng_only",label=TRUE, cols = big_region.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,order=levels(seuratObj),raster=FALSE)
dev.off()

##### class/subclass label transfer #####
# metainfo
zeng_anno <- readRDS("../../../03.data/02.metainfo/01.Young_Mouse/01-zeng_v3_metadata_downsample_1000.rds")

# subclass color
subclass.col = readRDS("../../../03.data/04.config_files/subclass_color.rds")

##### label transfer by cell type of max number in one seurat cluster #####
meta_zeng <- zeng_anno %>% column_to_rownames("cell_new")
rownames(seuratObj@meta.data[seuratObj@meta.data$group == "Zeng",]) <- gsub(".*@@_","",rownames(seuratObj@meta.data[seuratObj@meta.data$group == "Zeng",]))
meta_data.zeng <- seuratObj@meta.data[seuratObj@meta.data$group == "Zeng",]
meta_data.zeng$zeng_test <- rownames(meta_data.zeng)
meta_data.zeng$row <- rownames(meta_data.zeng)
meta_data.zeng$zeng_test <- gsub(".*@@_","",meta_data.zeng$zeng_test)
meta_zeng <- meta_zeng[rownames(meta_zeng) %in% meta_data.zeng$zeng_test,]
meta_zeng$row <- meta_data.zeng$row[match(rownames(meta_zeng),meta_data.zeng$zeng_test)]
rownames(meta_zeng) <- NULL
meta_zeng_final <- meta_zeng %>% column_to_rownames("row")
meta <- meta_zeng_final %>% dplyr::select(region_label,class_label,subclass_label)

Idents(seuratObj) <- seuratObj$seurat_clusters
levels(seuratObj)
class_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(seuratObj)) {
  cell_indices <- WhichCells(seuratObj, idents = ident_value)
  table_values <- table(meta[cell_indices, "class_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  class_df <- rbind(class_df, data.frame(idents = ident_value, max_variable = max_variable))
}

print(class_df)

seuratObj@meta.data$class_label = "NA"
for(i in 1:nrow(class_df)){
  seuratObj@meta.data[which(seuratObj@meta.data$seurat_clusters == class_df$idents[i]),'class_label'] <- class_df$max_variable[i]
}

## plot major class ##
## all ##
pdf("class_label.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label",label=FALSE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

pdf("class_label.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label",label=TRUE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

## Joint Cabernet ##
seuratObj@meta.data$class_label.color_Joint_Cabernet_only = seuratObj@meta.data$class_label
seuratObj@meta.data$class_label.color_Joint_Cabernet_only[which(seuratObj@meta.data$group == "Zeng")] = "Zeng"
Idents(seuratObj) <- seuratObj$class_label.color_Joint_Cabernet_only
levels(seuratObj)
pdf("class_label.color_Joint_Cabernet_only.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label.color_Joint_Cabernet_only",label=FALSE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(seuratObj)))
dev.off()

pdf("class_label.color_Joint_Cabernet_only.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label.color_Joint_Cabernet_only",label=TRUE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(seuratObj)))
dev.off()

## Zeng ##
seuratObj@meta.data$class_label.color_Zeng_only = seuratObj@meta.data$class_label
seuratObj@meta.data$class_label.color_Zeng_only[which(seuratObj@meta.data$group == "Joint_Cabernet")] = "Joint_Cabernet"
Idents(seuratObj) <- seuratObj$class_label.color_Zeng_only
levels(seuratObj)
pdf("class_label.color_Zeng_only.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label.color_Zeng_only",label=FALSE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = levels(seuratObj))
dev.off()

pdf("class_label.color_Zeng_only.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 15,height = 13)
DimPlot(seuratObj, group.by = "class_label.color_Zeng_only",label=TRUE, cols = major_class.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(seuratObj))
dev.off()

##### subclass #####
Idents(seuratObj) <- seuratObj$seurat_clusters
levels(seuratObj)
subclass_df <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(seuratObj)) {
  cell_indices <- WhichCells(seuratObj, idents = ident_value)
  table_values <- table(meta[cell_indices, "subclass_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  subclass_df <- rbind(subclass_df, data.frame(idents = ident_value, max_variable = max_variable))
}

print(subclass_df)

seuratObj@meta.data$subclass_label = "NA"
for(i in 1:nrow(subclass_df)){
  seuratObj@meta.data[which(seuratObj@meta.data$seurat_clusters == subclass_df$idents[i]),'subclass_label'] <- subclass_df$max_variable[i]
}

## plot subclass ##
## all ##
pdf("subclass_label.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label",label=FALSE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

pdf("subclass_label.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label",label=TRUE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE)
dev.off()

## Joint Cabernet ##
seuratObj@meta.data$subclass_label.color_Joint_Cabernet_only = seuratObj@meta.data$subclass_label
seuratObj@meta.data$subclass_label.color_Joint_Cabernet_only[which(seuratObj@meta.data$group == "Zeng")] = "Zeng"
Idents(seuratObj) <- seuratObj$subclass_label.color_Joint_Cabernet_only
levels(seuratObj)
pdf("subclass_label.color_Joint_Cabernet_only.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label.color_Joint_Cabernet_only",label=FALSE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(seuratObj)))
dev.off()

pdf("subclass_label.color_Joint_Cabernet_only.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label.color_Joint_Cabernet_only",label=TRUE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(seuratObj)))
dev.off()

## Zeng ##
seuratObj@meta.data$subclass_label.color_Zeng_only = seuratObj@meta.data$subclass_label
seuratObj@meta.data$subclass_label.color_Zeng_only[which(seuratObj@meta.data$group == "Joint_Cabernet")] = "Joint_Cabernet"
Idents(seuratObj) <- seuratObj$subclass_label.color_Zeng_only
levels(seuratObj)
pdf("subclass_label.color_Zeng_only.UMAP.without_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label.color_Zeng_only",label=FALSE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = levels(seuratObj))
dev.off()

pdf("subclass_label.color_Zeng_only.UMAP.with_label.label_transfer_by_cell_type_of_max_number_in_one_seurat_cluster.pdf",width = 20,height = 13)
DimPlot(seuratObj, group.by = "subclass_label.color_Zeng_only",label=TRUE, cols = subclass.col, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=levels(seuratObj))
dev.off()


saveRDS(seuratObj,file = "integration_Joint_Cabernet_and_Zeng.with_celltype.rds")

##### Joint Cabernet self cluster #####
Idents(seuratObj) = seuratObj@meta.data$group
Joint_Cabernet_all <- subset(seuratObj,cells = WhichCells(seuratObj,idents = "Joint_Cabernet"))
Joint_Cabernet <- readRDS("./rds/Joint_Cabernet_seurat.rds")
Joint_Cabernet <- FindClusters(Joint_Cabernet,resolution = 1.8,graph.name = 'SCT_snn')
Joint_Cabernet <- RunUMAP(Joint_Cabernet, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.3, n.neighbors=40, seed.use = 2)
pdf("seurat_clusters.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "seurat_clusters",label=FALSE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()

pdf("seurat_clusters.UMAP.with_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "seurat_clusters",label=TRUE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(Joint_Cabernet)))
dev.off()

meta_class_label <- Joint_Cabernet_all@meta.data %>% dplyr::select(class_label)
new_result_class_label <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(Joint_Cabernet)) {
  cell_indices <- WhichCells(Joint_Cabernet, idents = ident_value)
  table_values <- table(meta_class_label[cell_indices, "class_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  new_result_class_label <- rbind(new_result_class_label, data.frame(idents = ident_value, max_variable = max_variable))
}

print(new_result_class_label)
colnames(new_result_class_label) <- c("cluster","class_label")

Joint_Cabernet@meta.data$class_label = "NA"
for(i in 1:nrow(new_result_class_label)){
  Joint_Cabernet@meta.data[which(Joint_Cabernet@meta.data$seurat_clusters == new_result_class_label$cluster[i]),'class_label'] <- new_result_class_label$class_label[i]
}

pdf("class_label.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "class_label",label=FALSE, cols = major_class.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()

pdf("class_label.UMAP.with_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "class_label",label=TRUE, cols = major_class.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(Joint_Cabernet)))
dev.off()

meta_subclass_label <- Joint_Cabernet_all@meta.data %>% dplyr::select(subclass_label)
new_result_subclass_label <- data.frame(idents = character(), max_variable = character(), stringsAsFactors = FALSE)
for (ident_value in levels(Joint_Cabernet)) {
  cell_indices <- WhichCells(Joint_Cabernet, idents = ident_value)
  table_values <- table(meta_subclass_label[cell_indices, "subclass_label"])
  max_variable <- names(table_values)[which.max(table_values)]
  new_result_subclass_label <- rbind(new_result_subclass_label, data.frame(idents = ident_value, max_variable = max_variable))
}

print(new_result_subclass_label)
colnames(new_result_subclass_label) <- c("cluster","subclass_label")

Joint_Cabernet@meta.data$subclass_label = "NA"
for(i in 1:nrow(new_result_subclass_label)){
  Joint_Cabernet@meta.data[which(Joint_Cabernet@meta.data$seurat_clusters == new_result_subclass_label$cluster[i]),'subclass_label'] <- new_result_subclass_label$subclass_label[i]
}


pdf("subclass_label.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "subclass_label",label=FALSE, cols = subclass.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()

pdf("subclass_label.UMAP.with_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "subclass_label",label=TRUE, cols = subclass.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(Joint_Cabernet)))
dev.off()



##### label transfer once #####
Joint_Cabernet@meta.data$class_label_once = Joint_Cabernet_all@meta.data$class_label
Joint_Cabernet@meta.data$subclass_label_once = Joint_Cabernet_all@meta.data$subclass_label


pdf("class_label_once.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "class_label_once",label=FALSE, cols = major_class.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()
pdf("class_label_once.UMAP.with_label.Joint_Cabernet_self_cluster.pdf",width = 10,height = 7)
DimPlot(Joint_Cabernet, group.by = "class_label_once",label=TRUE, cols = major_class.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(Joint_Cabernet)))
dev.off()


pdf("subclass_label_once.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 17,height = 10)
DimPlot(Joint_Cabernet, group.by = "subclass_label_once",label=FALSE, cols = subclass.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()
pdf("subclass_label_once.UMAP.with_label.Joint_Cabernet_self_cluster.pdf",width = 17,height = 10)
DimPlot(Joint_Cabernet, group.by = "subclass_label_once",label=TRUE, cols = subclass.col,pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order=rev(levels(Joint_Cabernet)))
dev.off()


confusion_table_sub <- table(factor(Joint_Cabernet$class_label_once), factor(Joint_Cabernet$seurat_clusters))
confusion_matrix_normalize_sub <- apply(confusion_table_sub, 2, function(x) x / sum(x))
melted_confusion_matrix_normalize_sub <- melt(confusion_matrix_normalize_sub)
melted_confusion_matrix_normalize_sub = melted_confusion_matrix_normalize_sub[order(melted_confusion_matrix_normalize_sub$Var2,decreasing =T),]
melted_confusion_matrix_normalize_sub$Var2 = as.character(melted_confusion_matrix_normalize_sub$Var2)
melted_confusion_matrix_normalize_sub$Var2 = factor(melted_confusion_matrix_normalize_sub$Var2,levels = as.character(0:45))


pdf(paste0("confusion_heatmap.class_label_once.pdf"),height = 10,width=10)
ggplot(melted_confusion_matrix_normalize_sub, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "black", midpoint = 0.5, limits = c(0, 1)) + 
  theme_minimal(base_size = 20) + 
  labs(title="confusion matrix of Joint Cabernet RNA",x = "class_label_once", y = "seurat_cluster", fill = "Overlap Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 13), 
    axis.text.y = element_text(color = "black", size = 13), 
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"), 
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    legend.title = element_text(size = 15, face = "bold"), 
    legend.text = element_text(size = 15), 
    panel.border = element_rect(color = "black", fill = NA, size = 1.5), 
    axis.ticks.length = unit(0.3, "cm"), 
    axis.ticks = element_line(color = "black", linewidth = 1.5) 
  )
dev.off()

confusion_table_sub <- table(factor(Joint_Cabernet$subclass_label_once), factor(Joint_Cabernet$seurat_clusters))
confusion_matrix_normalize_sub <- apply(confusion_table_sub, 2, function(x) x / sum(x))
melted_confusion_matrix_normalize_sub <- melt(confusion_matrix_normalize_sub)
melted_confusion_matrix_normalize_sub = melted_confusion_matrix_normalize_sub[order(melted_confusion_matrix_normalize_sub$Var2,decreasing =T),]
melted_confusion_matrix_normalize_sub$Var2 = as.character(melted_confusion_matrix_normalize_sub$Var2)
melted_confusion_matrix_normalize_sub$Var2 = factor(melted_confusion_matrix_normalize_sub$Var2,levels = as.character(0:45))


pdf(paste0("confusion_heatmap.subclass_label_once.pdf"),height = 10,width=15)
ggplot(melted_confusion_matrix_normalize_sub, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "red", high = "black", midpoint = 0.5, limits = c(0, 1)) +  
  theme_minimal(base_size = 20) +  
  labs(title="confusion matrix of Joint Cabernet RNA",x = "subclass_label_once", y = "seurat_cluster", fill = "Overlap Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 13), 
    axis.text.y = element_text(color = "black", size = 13),  
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),  
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  
    legend.title = element_text(size = 15, face = "bold"),  
    legend.text = element_text(size = 15),  
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),  
    axis.ticks.length = unit(0.3, "cm"),  
    axis.ticks = element_line(color = "black", linewidth = 1.5)  
  )
dev.off()


saveRDS(Joint_Cabernet,file = "Joint_Cabernet.with_celltype.rds")

##### 
confusion_table_sub <- table(factor(Joint_Cabernet$subclass_label_once), factor(Joint_Cabernet$seurat_clusters))
confusion_matrix_normalize_sub <- apply(confusion_table_sub, 2, function(x) x / sum(x))

pdf("subclass_confusion_matrix_of_label_transfer_by_max_abundant_cell_type.heatmap.pdf",width=15,height=10)
col_fun =colorRamp2(breaks=c(0,0.5,1),color=c("white", "red", "black"))
Heatmap(t(confusion_matrix_normalize_sub),
        column_gap = unit(0, "points"),
        col = col_fun,
        cluster_rows=T,
        cluster_columns=F,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names=T,
        show_row_names=T,
        row_names_side="left",
        column_names_side="bottom",
        column_names_rot = 75,
        column_dend_side = "top",
        border = TRUE,
        row_title = "subclass",
        column_title = "Seurat cluster",
        column_title_gp = gpar(fontsize = 25, fontface = "bold"),
        column_title_side = "top",
        column_title_rot = 0,
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 25, fontface = "bold"),
        heatmap_legend_param = list(title_position = "lefttop-rot",title = "Overlap Score",legend_height = unit(4, "cm"),fontsize = 20),
        use_raster = FALSE)
dev.off()



confusion_table_sub <- table(factor(Joint_Cabernet$subclass_label_once), factor(Joint_Cabernet$seurat_clusters))
confusion_matrix_normalize_sub <- apply(confusion_table_sub, 2, function(x) x / sum(x))
write.csv(data.frame(seurat_cluster = names(apply(confusion_matrix_normalize_sub, 2, function(x) max(x))), max_ratio = apply(confusion_matrix_normalize_sub, 2, function(x) max(x))),file = "seurat_cluster_max_ratio_of_label_transfer_by_max_abundant_cell_type.csv",quote=F,row.names=F,col.names=T)


test = Joint_Cabernet
test <- FindClusters(test)
test <- RunUMAP(test, dims=1:30, dim.embed=5, reduction="pca", min.dist=0.4, n.neighbors=80)

pdf("test_seurat.UMAP.without_label.Joint_Cabernet_self_cluster.pdf",width = 13,height = 7)
DimPlot(test, group.by = "seurat_clusters",label=FALSE, pt.size =0.5,seed=1100,label.box=T,label.size = 4,raster=FALSE,order = rev(levels(Joint_Cabernet)))
dev.off()




