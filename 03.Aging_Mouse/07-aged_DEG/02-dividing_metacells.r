#########    All "our" in the following code refers to Joint Cabernet.
library(hdWGCNA,lib.loc="/share/home/liuyy/anaconda3/envs/liuyy_R/lib/R/library")

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

#####The results of hdWGCNA are saved in the Seurat object@misc slot
our_seurat<-readRDS(paste0(indir,"/our_seruat_with_cluster_corrected.rds"))
DefaultAssay(our_seurat) <- "RNA"

##Use all genes after count filtering
filtered_genes<-read.csv("../../output/03-aging/07-aged_DEG/count_cpm_filtered_gene.csv")
our_seurat <- SetupForWGCNA(
  our_seurat,
  gene_select = "custom", # the gene selection approach
  features=filtered_genes$x,
  wgcna_name = "tutorial" #Specify the hdWGCNA experiment name
)
length(our_seurat@misc$tutorial$wgcna_genes)


######## dividing metacells
#Divide metacells separately for each subclass and age

our_seurat <- MetacellsByGroups(
  seurat_obj = our_seurat,
  group.by = c("lt_twice_subclass", "age"), # specify the columns in seurat_obj@meta.data to group by
  k = 15, # Number of nearest neighbors
  max_shared = 3, # Both Metacells share maximum cells
  min_cells=40,  #Minimum number of cells per metacell, default is 100.
  ident.group = 'lt_twice_subclass' # set the Idents of the metacell seurat object.ident.group must be in group.by.
)

metacell_obj <- GetMetacellObject(our_seurat)
metacell_obj$subclass_age<-paste0(metacell_obj$lt_twice_subclass,"-",metacell_obj$age)
metacell_obj <- NormalizeData(metacell_obj)
saveRDS(metacell_obj,paste0(outdir,"/metacells_obj_k15_mincell40_maxsh3_corrected.rds"))





