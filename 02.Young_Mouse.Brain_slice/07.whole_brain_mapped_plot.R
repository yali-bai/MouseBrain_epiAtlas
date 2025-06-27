#########    All "Joint_cabernet" in the following code refers to Joint Cabernet.
##### 01.import packages #####
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

##### 01.extract loci of background #####
spatial_zhuang_merfish_loci.df = read.csv("../04.data/03.download_data/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
spatial_zhuang_merfish_loci.df = spatial_zhuang_merfish_loci.df[intersect(which(spatial_zhuang_merfish_loci.df$z > 7.33),which(spatial_zhuang_merfish_loci.df$z < 7.34)),] # subset cells by Z axis around 7.33 which is consistent with Joint-Cabernet MERFISH data
rownames(spatial_zhuang_merfish_loci.df) = spatial_zhuang_merfish_loci.df$cell_label
spatial_zhuang_merfish_loci.df = spatial_zhuang_merfish_loci.df[,c("x","y")]
colnames(spatial_zhuang_merfish_loci.df) = c("umap_1", "umap_2")

##### 02.color Exc neuron cells on whole brain #####
merged.seuratobj.sct = readRDS("../output/02.Young_Mouse.Brain_slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
merged.seuratobj.sct_subset = subset(merged.seuratobj.sct,three_class == "Exc") # subset by three class information 
plot_df = data.frame(Embeddings(merged.seuratobj.sct_subset, reduction = 'umap')[rownames(merged.seuratobj.sct_subset@meta.data)[which(merged.seuratobj.sct_subset$source == "Joint_cabernet")],])
plot_df$subclass = merged.seuratobj.sct@meta.data[rownames(plot_df),"subclass"]
inte.col <- c("Joint_cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#89C75F", "Astro-TE NN" = "#0C727C",
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

pdf("../output/02.Young_Mouse.Brain_slice/map/Exc_spatial.whole_biran.pdf",width = 40,height = 50)
ggplot() +
            geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") + # size = 15
            geom_point(data = plot_df, mapping = aes(umap_1, umap_2, color=subclass), size = 12) +
            scale_color_manual(values=inte.col)+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) +
            labs(title = "subclass", x = "spatial_1", y = "spatial_2")
dev.off()

##### 03.color Inh neuron cells on whole brain #####
merged.seuratobj.sct = readRDS("../output/02.Young_Mouse.Brain_slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
merged.seuratobj.sct_subset = subset(merged.seuratobj.sct,three_class == "Inh")
plot_df = data.frame(Embeddings(merged.seuratobj.sct_subset, reduction = 'umap')[rownames(merged.seuratobj.sct_subset@meta.data)[which(merged.seuratobj.sct_subset$source == "Joint_cabernet")],])
plot_df$subclass = merged.seuratobj.sct@meta.data[rownames(plot_df),"subclass"]
inte.col <- c("Joint_cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#89C75F", "Astro-TE NN" = "#0C727C",
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

pdf("../output/02.Young_Mouse.Brain_slice/map/Inh_spatial.whole_biran.pdf",width = 40,height = 50)
ggplot() +
            geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") + # size = 15
            geom_point(data = plot_df, mapping = aes(umap_1, umap_2, color=subclass), size = 12) +
            scale_color_manual(values=inte.col)+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) +
            labs(title = "subclass", x = "spatial_1", y = "spatial_2")
dev.off()

##### 04.color Non-neuron cells on whole brain #####
merged.seuratobj.sct = readRDS("../output/02.Young_Mouse.Brain_slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
merged.seuratobj.sct_subset = subset(merged.seuratobj.sct,three_class == "Non-neuron")
plot_df = data.frame(Embeddings(merged.seuratobj.sct_subset, reduction = 'umap')[rownames(merged.seuratobj.sct_subset@meta.data)[which(merged.seuratobj.sct_subset$source == "Joint_cabernet")],])
plot_df$subclass = merged.seuratobj.sct@meta.data[rownames(plot_df),"subclass"]

inte.col <- c("Joint_cabernet"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = "#89C75F", "Astro-TE NN" = "#0C727C",
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
                               
pdf("../output/02.Young_Mouse.Brain_slice/map/Non_neuron_spatial.whole_biran.pdf",width = 40,height = 50)
ggplot() +
            geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") + # size = 15
            geom_point(data = plot_df, mapping = aes(umap_1, umap_2, color=subclass), size = 12) +
            scale_color_manual(values=inte.col)+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) +
            labs(title = "subclass", x = "spatial_1", y = "spatial_2")
dev.off()