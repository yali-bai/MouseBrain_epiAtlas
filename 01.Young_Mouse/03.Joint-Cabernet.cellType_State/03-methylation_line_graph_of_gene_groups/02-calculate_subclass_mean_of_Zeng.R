library(Seurat)

seuratObj = readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/integration_Joint_Cabernet_and_Zeng.with_celltype.rds")

Zeng.obj = subset(seuratObj,group == "Zeng")

metainfo = seuratObj@meta.data

for(i in 1:length(unique(Zeng.obj@meta.data$subclass_label))){
    print(paste0(i,": ",unique(Zeng.obj@meta.data$subclass_label)[i]))
    if(i ==1 ){
        result.df = data.frame(temp=rowMeans(Zeng.obj@assays$RNA$data[,Zeng.obj@meta.data$subclass_label == unique(Zeng.obj@meta.data$subclass_label)[i]]))
        colnames(result.df) = unique(Zeng.obj@meta.data$subclass_label)[i]
    }
    else{
        result.df = cbind(result.df,data.frame(temp = rowMeans(Zeng.obj@assays$RNA$data[,Zeng.obj@meta.data$subclass_label == unique(Zeng.obj@meta.data$subclass_label)[i]])))
        colnames(result.df)[dim(result.df)[2]] = unique(Zeng.obj@meta.data$subclass_label)[i]
    }
}

write.csv(result.df,file = "zeng_subclass_mean_dat_final.csv",quote=F,row.names=T,col.names=T)
