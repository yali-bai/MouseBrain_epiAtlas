library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(gridExtra)
library(ggplot2)

seuratObj <- readRDS("../../01.RNA-integration/04.Joint-Cabernet.Zeng_10X_RNA.integration/Joint_Cabernet.with_celltype.rds")

paired_sampleinfo = read.csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",header =T)

DefaultAssay(seuratObj) = "RNA"
seuratObj = NormalizeData(seuratObj)

RNA_mat = seuratObj@assays$RNA$data %>% data.frame()
colnames(RNA_mat) = unlist(lapply(colnames(RNA_mat), function(x) strsplit(x,'Cabernet.._')[[1]][2]))
rownames(RNA_mat) = unlist(lapply(rownames(RNA_mat), function(x) strsplit(x,'\\.')[[1]][1]))

RNA_mat = RNA_mat[,paired_sampleinfo[paired_sampleinfo$total_QC == 1,"RNA_SampleID"]]
colnames(RNA_mat) = paired_sampleinfo[paired_sampleinfo$total_QC == 1,"unique_id"]

RNA_mat <- t(RNA_mat)

UMAP <- data.frame(Embeddings(seuratObj, reduction = 'umap')[Cells(seuratObj),])
rownames(UMAP) = unlist(lapply(rownames(UMAP), function(x) strsplit(x,'Cabernet@@_')[[1]][2]))
UMAP = UMAP[paired_sampleinfo[paired_sampleinfo$total_QC == 1,"RNA_SampleID"],]
rownames(UMAP) = paired_sampleinfo[paired_sampleinfo$total_QC == 1,"unique_id"]


mm10_genes <- read.delim('../../../03.data/01.ref/mm10.genes.bed',header = T)
rownames(mm10_genes) = unlist(lapply(mm10_genes$gene_id, function(x) strsplit(x,'\\.')[[1]][1]))
colnames(RNA_mat) <- mm10_genes[colnames(RNA_mat),"gene_name"]

RNA_final <- cbind(UMAP,RNA_mat)


RNA_plot = ggplot(RNA_final[,c("umap_1","umap_2","Cdh18")], aes(x = umap_1, y = umap_2, color = Cdh18)) +
    geom_point(size=0.1) +
    labs(x = "", y = "",title = "RNA") +
    theme_minimal() +
    theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
          axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
    scale_color_gradient(low = "lightgrey", high = "#89288F")

read_mat <- function(group,type,mc_type){
    if(type == "5hmC"){
        coln = "hmC_SampleID"
    }else{
        coln = "mC_SampleID"
    }
    data <- readRDS(sprintf("temp/%s_geneslop2k.%s.fill_na.rds",group,mc_type))
    data <- t(data)
    colnames(data) <- mm10_genes[colnames(data),"gene_name"]
    rownames(data) = str_replace_all(rownames(data),"allc_","")
    data = data[paired_sampleinfo[paired_sampleinfo$total_QC == 1,coln],]
    data <- cbind(UMAP,data)
    colnames(data) = c("umap_1","umap_2","Cdh18")
    return(data)
}

hmCpG_final <- read_mat("5hmC","5hmC","CG")
mCpG_final <- read_mat("5mC_5hmC","5mC","CG")
true_CpG_final <- read_mat("5mC","5mC","CG")
true_CH_final <- read_mat("5mC","5mC","CH")

plot_gene_expression <- function(gene_name,RNA_low,RNA_high) {
    
    data_frames <- list(
        RNA_final = RNA_final,
        hmCpG_final = hmCpG_final,
        mCpG_final = mCpG_final,
        true_CpG_final = true_CpG_final,
        true_CH_final = true_CH_final
    )

    for (df_name in names(data_frames)) {
        if (!(gene_name %in% colnames(data_frames[[df_name]]))) {
            warning(paste("Gene", gene_name, "not found in", df_name))
            return(NULL)
        }
    }
    
    RNA_final[[gene_name]] <- pmin(pmax(RNA_final[[gene_name]], RNA_low), RNA_high)
    hmCpG_final[[gene_name]] <- pmin(pmax(hmCpG_final[[gene_name]], 0.25), 0.40)
    mCpG_final[[gene_name]] <- pmin(pmax(mCpG_final[[gene_name]], 0.3), 1)
    true_CpG_final[[gene_name]] <- pmin(pmax(true_CpG_final[[gene_name]], 0.3), 1)
    true_CH_final[[gene_name]] <- pmin(pmax(true_CH_final[[gene_name]], 0), 0.05)
    
    p1 <- ggplot(RNA_final[,c("umap_1","umap_2",gene_name)], aes_string(x = "umap_1", y = "umap_2", color = gene_name)) +
        geom_point(size=0.3) +
        labs(x = "", y = "",title = "RNA") +
        theme_minimal() +
        theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
              axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
        scale_color_gradient(low = "lightgrey", high = "#89288F",limits = c(RNA_low, RNA_high),breaks = c(RNA_low, RNA_high))

    pp1 <- p1 + theme(legend.position = "none")
    legend_p1 <- get_legend(p1)

    p2 <- ggplot(hmCpG_final[,c("umap_1","umap_2",gene_name)], aes_string(x = "umap_1", y = "umap_2", color = gene_name)) +
        geom_point(size=0.3) +
        labs(x = "", y = "",title = "5hmCG") +
        theme_minimal() +
        theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
              axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
        scale_color_gradient(low = "darkblue", high = "#e6c200",limits = c(0.25, 0.40),breaks = c(0.25, 0.40))

    pp2 <- p2 + theme(legend.position = "none")
    legend_p2 <- get_legend(p2)

    p3 <- ggplot(mCpG_final[,c("umap_1","umap_2",gene_name)], aes_string(x = "umap_1", y = "umap_2", color = gene_name)) +
        geom_point(size=0.3) +
        labs(x = "", y = "",title = "5mCG+5hmCG") +
        theme_minimal() +
        theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
              axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
        scale_color_gradient(low = "darkblue", high = "#e6c200",limits = c(0.3, 1),breaks = c(0.3, 1))

    pp3 <- p3 + theme(legend.position = "none")
    legend_p3 <- get_legend(p3)

    p4 <- ggplot(true_CpG_final[,c("umap_1","umap_2",gene_name)], aes_string(x = "umap_1", y = "umap_2", color = gene_name)) +
        geom_point(size=0.3) +
        labs(x = "", y = "",title = "5mCG") +
        theme_minimal() +
        theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
              axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
        scale_color_gradient(low = "darkblue", high = "#e6c200",limits = c(0.3, 1),breaks = c(0.3, 1))

    pp4 <- p4 + theme(legend.position = "none")
    legend_p4 <- get_legend(p4)

    p5 <- ggplot(true_CH_final[,c("umap_1","umap_2",gene_name)], aes_string(x = "umap_1", y = "umap_2", color = gene_name)) +
        geom_point(size=0.3) +
        labs(x = "", y = "",title = "5mCH") +
        theme_minimal() +
        theme(text = element_text(size = 15),panel.grid = element_blank(),axis.text = element_blank(),
              axis.ticks = element_blank(),plot.title = element_text(vjust = 0.5,hjust = 0.5,size=40)) +
        scale_color_gradient(low = "darkblue", high = "#75D054FF",limits = c(0, 0.05),breaks = c(0, 0.05))

    pp5 <- p5 + theme(legend.position = "none")
    legend_p5 <- get_legend(p5)
    
    gpp1 <- arrangeGrob(pp1, heights = unit(6.2, "inches"), widths = unit(6, "inches"))
    gpp2 <- arrangeGrob(pp2, heights = unit(6.2, "inches"), widths = unit(6, "inches"))
    gpp3 <- arrangeGrob(pp3, heights = unit(6.2, "inches"), widths = unit(6, "inches"))
    gpp4 <- arrangeGrob(pp4, heights = unit(6.2, "inches"), widths = unit(6, "inches"))
    gpp5 <- arrangeGrob(pp5, heights = unit(6.2, "inches"), widths = unit(6, "inches"))
    
    glegend_p1 <- arrangeGrob(legend_p1, heights = unit(3, "inches"), widths = unit(3, "inches"))
    glegend_p2 <- arrangeGrob(legend_p2, heights = unit(3, "inches"), widths = unit(3, "inches"))
    glegend_p3 <- arrangeGrob(legend_p3, heights = unit(3, "inches"), widths = unit(3, "inches"))
    glegend_p4 <- arrangeGrob(legend_p4, heights = unit(3, "inches"), widths = unit(3, "inches"))
    glegend_p5 <- arrangeGrob(legend_p7, heights = unit(3, "inches"), widths = unit(3, "inches"))

    combined_plot <- arrangeGrob(gpp1, gpp3, gpp2, gpp4, gpp5, 
                                glegend_p1, glegend_p3, glegend_p2, glegend_p4, glegend_p5, ncol = 5, nrow = 2)

    return(combined_plot)
}



ggsave("Cdh18.pdf", 
       plot_gene_expression("Cdh18",0,5), 
       width = 32, height = 13)