#########    All "our" in the following code refers to Joint Cabernet.
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
library(data.table) 

##### 02.change working path #####
# setwd("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/MERFISH/20240902")

##### 03.load data after filling na #####
#load("RNA_DNA_fill_na.20240925.RData")
load("../output/02-slice/RNA_DNA_fill_na.20240925.RData")

##### 04.data prepare #####
## subclass order ##
#subclass_order.df = read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/sample_info/04_order_for_class_subclass/subclass_order_for_integration_with_zhuang.txt",header=F)
subclass_order.df = read.csv("../input/02-slice/subclass_order_for_integration_with_zhuang.txt",header=F)
order.v = subclass_order.df$V1
## paired information of RNA and DNA sample id ##
#paired_sampleinfo = read.csv("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/sample_info/01_Sampleinfo/RNA_DNA_match_name_QC_MERFISH.csv",header=T)
paired_sampleinfo = read.csv("../input/02-slice/RNA_DNA_match_name_QC_MERFISH.csv",header=T)
## loci and subclass information ##
#merged.seuratobj.sct <- readRDS("merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
merged.seuratobj.sct <- readRDS("../output/02-slice/map/merged.seuratobj.sct.loci_transfer.the_nearst_1_cell.rds")
metainfo = merged.seuratobj.sct@meta.data
## Zhuang MERFISH metainfo ##
#loc.df <- read.csv("/share/analysisdata/Methyl/public/analysis/data/MERFISH/Zhuang_dataset/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
loc.df <- read.csv("../input/02-slice/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
rownames(loc.df) <- loc.df$cell_label
TSNE_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'tsne')[Cells(merged.seuratobj.sct),])  # extract tsne loci for substitution
spatial_zhuang_merfish_loci.df <- as.data.frame(TSNE_RNA[which(is.na(str_match(rownames(TSNE_RNA),"Mouses"))),])
spatial_zhuang_merfish_loci.df$spatial_1 <- loc.df[rownames(spatial_zhuang_merfish_loci.df),"x"]
spatial_zhuang_merfish_loci.df$spatial_2 <- loc.df[rownames(spatial_zhuang_merfish_loci.df),"y"]
## loci after mapping ##
TSNE_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  # extract final spatial loci

## color ##
## for DNA ##
colors <- colorRampPalette(c("#253494", "#225EA8", "#1D91C0", "#41B6C4","#FFED6F","#FFFF33","#FFFF00")) # blue to yellow
color_gradient <- colors(75)

## for RNA ##
colors <- colorRampPalette(c("#5614b0", "#7438bd","#915cca","#af80d6","#cca4e3","#F4C800","#FFED6F","#FFFF00")) ## purple to yellow
color_gradient_RNA <- colors(75)

##### 05. plot #####
## define plot function ##
## 01. plot DNA methyl levels of one gene spatially ##
marker_plot_limits <- function(cell_name,methyl_matrix,name,gene,min,max){
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
        rownames(paired_sampleinfo) = paired_sampleinfo$RNA
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
                limits = c(min, max), # limit legend range
                breaks = seq(min,max,(max-min)/2),
                colors = color_gradient,
                na.value = color_gradient[length(color_gradient)])+
            theme_minimal() +
            theme(text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) +
            labs(title = name, x = "spatial_1", y = "spatial_2") 
}

## rotate spatial loci 90 degrees clockwise ##
spatial_loci <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  # extract final tsne loci
spatial_loci = as.matrix(spatial_loci*-1)
merged.seuratobj.sct@reductions$umap@cell.embeddings <- spatial_loci
TSNE_RNA = TSNE_RNA*-1
spatial_zhuang_merfish_loci.df = spatial_zhuang_merfish_loci.df*-1
spatial_zhuang_merfish_loci.df$umap_1 = spatial_zhuang_merfish_loci.df$spatial_1
spatial_zhuang_merfish_loci.df$umap_2 = spatial_zhuang_merfish_loci.df$spatial_2

## 02. define plot function which plot celltype, Joint-Cabernet RNA, Zhuang MERFISH RNA, Joint-Cabernet DNA ##
limits_marker_RNA_cluster_plot_neuron_type <- function(marker_type, type, feature, cl, gene, pos_neg,max_vector){
                if(cl == "L2/3 IT CTX Glut"){
                    ## define color ##
                    inte.col <- c("Our"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = 'lightgrey', "Astro-TE NN" = 'lightgrey',
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
                    ## plot subclass ##
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
                    ## define color ##
                    inte.col <- c("Our"='lightgrey',"Zhuang"='lightgrey',"Oligo NN" = 'lightgrey', "Astro-TE NN" = 'lightgrey',
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
                        ## plot subclass ##
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

                ## extract RNA info: expression level and loci ##
                plot_df <- as.data.frame(cbind(as.numeric(RNA_fill[strsplit(gene,"\\.")[[1]][1],which(!is.na(str_match(colnames(RNA_fill),"Mouse")))]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))],"umap_1"]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))],"umap_2"])))
                colnames(plot_df) <- c("RNA","umap_1","umap_2")
                rownames(plot_df) <- colnames(merged.seuratobj.sct@assays$RNA$data)[which(!is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))]

                plot_df <- plot_df[cell_n.v,]
                plot_df <- na.omit(plot_df)
                plot_df[which(plot_df$RNA < max_vector[1]),"RNA"] = max_vector[1]
                ## background loci ##
                spatial_zhuang_merfish_loci.df <- TSNE_RNA[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse"))),]
                ## Joint-Cabernet RNA plot ##
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
                labs(title = "RNA.Our", x = "spatial_1", y = "spatial_2")

                ## extract Zhuang MERFISH cells of the same subclasses ##
                zhuang_plot_df <- as.data.frame(cbind(as.numeric(RNA_fill[strsplit(gene,"\\.")[[1]][1],which(is.na(str_match(colnames(RNA_fill),"Mouse")))]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))],"umap_1"]),
                    as.numeric(TSNE_RNA[colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))],"umap_2"])))
                colnames(zhuang_plot_df) <- c("RNA","umap_1","umap_2")
                rownames(zhuang_plot_df) <- colnames(merged.seuratobj.sct@assays$RNA$data)[which(is.na(str_match(colnames(merged.seuratobj.sct@assays$RNA$data),"Mouse")))]
                zhuang_plot_df <- zhuang_plot_df[cell_n.v,]
                zhuang_plot_df <- na.omit(zhuang_plot_df)
                zhuang_plot_df[which(zhuang_plot_df$RNA < max_vector[3]),"RNA"] = max_vector[3]
                ## Zhuang MERFISH RNA plot ##
                zhuang_RNA_plot <- ggplot() + 
                geom_point(data = spatial_zhuang_merfish_loci.df, mapping = aes(umap_1, umap_2), size = 3, color = "lightgrey") +
                geom_point(data = zhuang_plot_df, mapping = aes(umap_1, umap_2, color=RNA), size = 6) +
                scale_color_gradientn(  
                    limits = c(max_vector[3], max_vector[4]), # limit legend range 
                    breaks = seq(max_vector[3], max_vector[4],1),
                    colors = color_gradient_RNA,
                    na.value = color_gradient_RNA[length(color_gradient_RNA)])+
                theme_minimal() +
                theme(text = element_text(face="bold", size=30),panel.grid = element_blank(),axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_text(face="bold", size=30),
                    axis.title.y = element_text(face="bold", size=30))+
                labs(title = "RNA.Zhuang", x = "spatial_1", y = "spatial_2") 

                ## define DNA ##
                if(feature == "genebody"){
                    ## mCG
                    mC_CGN = geneslop2k_mC_CGN_fill
                                    
                    ## hmCG
                    hmC_CGN = geneslop2k_hmC_CGN_fill
                    
                    ## mCG_hmCG
                    hmCG_mCG = geneslop2k_hmCG_mCG_fill
                    
                    ## mCH
                    mC_CHN = geneslop2k_mC_CHN_fill
                    
                    ## hmCH
                    hmC_CHN = geneslop2k_hmC_CHN_fill
                    
                    ## mCH_hmCH
                    hmCH_mCH = geneslop2k_hmCH_mCH_fill
                }else if(feature == "promoter"){
                    ## mCG
                    mC_CGN = promoter_mC_CGN_fill
                
                    ## hmCG
                    hmC_CGN = promoter_hmC_CGN_fill

                    ## mCG_hmCG
                    hmCG_mCG = promoter_hmCG_mCG_fill

                    ## mCH
                    mC_CHN = promoter_mC_CHN_fill

                    ## hmCH
                    hmC_CHN = promoter_hmC_CHN_fill

                    ## mCH_hmCH
                    hmCH_mCH = promoter_hmCH_mCH_fill
                }

                for(i in 1:dim(mC_CGN)[2]){
                    mC_CGN[which(mC_CGN[,i] < 0),i] <- 0
                }
                for(i in 1:dim(mC_CHN)[2]){
                    mC_CHN[which(mC_CHN[,i] < 0),i] <- 0
                }

                ## run plot function ##
                mCG_plot <- marker_plot_limits(cell_n.v,mC_CGN,"5mCG",strsplit(gene,"\\.")[[1]][1],max_vector[5],max_vector[6])
                hmCG_plot <- marker_plot_limits(cell_n.v,hmC_CGN,"5hmCG",strsplit(gene,"\\.")[[1]][1],max_vector[7],max_vector[8])
                mCG_hmCG_plot <- marker_plot_limits(cell_n.v,hmCG_mCG,"5mCG+5hmCG",strsplit(gene,"\\.")[[1]][1],max_vector[9],max_vector[10])
                mCH_plot <- marker_plot_limits(cell_n.v,mC_CHN,"5mCH",strsplit(gene,"\\.")[[1]][1],max_vector[11],max_vector[12])
                hmCH_plot <- marker_plot_limits(cell_n.v,hmC_CHN,"5hmCH",strsplit(gene,"\\.")[[1]][1],max_vector[13],max_vector[14])
                mCH_hmCH_plot <- marker_plot_limits(cell_n.v,hmCH_mCH,"5mCH+5hmCH",strsplit(gene,"\\.")[[1]][1],max_vector[15],max_vector[16])

                ## arrange the figures ##
                plot_grid(  
                    plotlist = list(  
                        plot_grid(all_subclass_plot, RNA_plot, zhuang_RNA_plot, ncol = 3, rel_widths = c(8,8,8)),  
                        plot_grid(mCG_plot, hmCG_plot, mCG_hmCG_plot, ncol = 3, rel_widths = c(8,8,8)),  
                        plot_grid(mCH_plot, hmCH_plot, mCH_hmCH_plot, ncol = 3, rel_widths = c(8,8,8))
                    ),  
                    nrow = 3
                )

}

## run function ##
plot <- limits_marker_RNA_cluster_plot_neuron_type("RNA","subclass","genebody","L2/3 IT CTX Glut","ENSMUSG00000001985.9","positive",c(0,2,0.5,4.5,0.35,0.65,0.25,0.4,0.5,0.9,0.01,0.04,0.005,0.01,0.025,0.04))
pdf("../output/02-slice/pseudotime/IT_Glut.ENSMUSG00000001985.9.genebody.pdf",width = 50,height = 40)
print(plot)
dev.off()