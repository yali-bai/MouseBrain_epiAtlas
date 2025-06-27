##### 01. import packages #####
library(Seurat)
library(ggplot2)
library(dplyr)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 03. read files in #####
our.seuratobj = readRDS(paste0(indir,"our.seuratobj.min_dist_0.5.rds"))

global_hmCG_diff = read.csv("../../04.data/05.intermediate_files/02.DNA/02.Aging_Mouse/global_DNA_old_vs_young_diff.CG.csv",header=T)
rownames(global_hmCG_diff) = global_hmCG_diff$lt_twice_subclass
our.seuratobj$global_hmCG_diff = global_hmCG_diff[our.seuratobj@meta.data$lt_twice_subclass,"hmCG_diff"]

umap_data <- data.frame(Embeddings(our.seuratobj, reduction = 'umap')[Cells(our.seuratobj),])
umap_data$global_5hmCG_diff = (our.seuratobj$global_hmCG_diff)


pdf(paste0("../../output/03.Aging_Mouse/global_5hmCG_diff_in_subclass.umap.without_legend.pdf"),width = 47.569,height = 43.592)
ggplot() +
            geom_point(data = umap_data, mapping = aes(umap_1, umap_2, color = global_5hmCG_diff), size = 15) + # size = 15 
            scale_color_gradientn(colors = c("#b7b7b9","#c1c1c3","#D3D3D3","#746ea3","#61599d"),breaks = c(0.01,0.015,0.02,0.025,0.03))+ # grey to purple
            theme_minimal() +
            theme(legend.position="none",
                #legend.key.size = unit(45, "pt"),
                #legend.title = element_text(face="bold",size=30),
                #legend.text = element_text(face="bold",size=30),
                text = element_text(face="bold",size = 30),
                panel.grid = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_text(face="bold", size=30),
                axis.title.y = element_text(face="bold", size=30)) #+
            #labs(title = "UMAP Plot with Coloring by subclass mean global genebody 5hmCG diff", x = "UMAP 1", y = "UMAP 2")
dev.off()






