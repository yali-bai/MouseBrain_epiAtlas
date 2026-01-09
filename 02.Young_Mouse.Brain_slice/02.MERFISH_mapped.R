library(getopt)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

arguments = matrix(c(
  'help', 'h', 0, "logical", "",
  'neuron_type', 'n', 1, "character", "",
  'label', 'l', 1, "character", ""
), byrow=TRUE, ncol=5)
args = getopt(arguments)

# if help was asked for print a friendly message
# and exit with a non-zero error code
#if ( !is.null(args$help) ) {
#  cat(getopt(arguments, usage=TRUE))
#  q(status=1)
#}

if (!is.null(args$help) || is.null(args$neuron_type) || is.null(args$label)) {
  cat(paste(getopt(arguments, usage = T), "\n"))
  q()
}


merged.seuratobj.sct <- readRDS(paste0(args$neuron_type,"_obj.label_transfer_",args$label,".rds"))

##### extract umap dis
UMAP_RNA <- data.frame(Embeddings(merged.seuratobj.sct, reduction = 'umap')[Cells(merged.seuratobj.sct),])  # extract final tsne loci

##### calculate 
dist_matrix <- as.dist(dist(UMAP_RNA))
dist_matrix_mat <- as.matrix(dist_matrix)  

dist_matrix_mat <- dist_matrix_mat[which(!is.na(str_match(rownames(dist_matrix_mat),"TSO"))),]   # reserve our RNA sample as row
dist_matrix_mat <- dist_matrix_mat[,which(is.na(str_match(colnames(dist_matrix_mat),"TSO")))]   # reserve MERFISH sample as col

metainfo <- merged.seuratobj.sct@meta.data
metainfo_subset <- metainfo[which(is.na(str_match(rownames(metainfo),"TSO"))),]

##### calculate nearst MERFISH cell
result_df <- data.frame(idents = character(), nearst = character(), stringsAsFactors = FALSE)
for(i in 1:dim(dist_matrix_mat)[1]){
    cl <- metainfo[rownames(dist_matrix_mat)[i],"subclass"]
    region <- metainfo[rownames(dist_matrix_mat)[i],"major_brain_region_v2"]
    cells <- rownames(metainfo[intersect(which(metainfo$subclass == cl),which(metainfo$major_brain_region_v2 == region)),])
    temp <- dist_matrix_mat[i,match(cells,colnames(dist_matrix_mat))]
    indices <- order(temp,decreasing = FALSE)  # extract the nearst 10 sample index
    result_df <- rbind(result_df, data.frame(idents = rownames(dist_matrix_mat)[i], nearst = paste(names(temp)[indices[1:1]],collapse=";"))) # record
}

# Zhuang MERFISH metainfo including cell_label, loci 
loc.df <- read.csv("../03.data/03.download_data/cell_metadata_Zhuang_MERFISH.csv",colClasses = c("character","character","character","character","character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
rownames(loc.df) <- loc.df$cell_label

x.v <- c()
y.v <- c()
for(i in 1:dim(result_df)[1]){
    cells <- unlist(lapply(result_df[i,"nearst"], function(x) strsplit(x,";")[[1]]))
    x.v <- c(x.v,mean(loc.df[cells,"x"]))
    y.v <- c(y.v,mean(loc.df[cells,"y"]))   # calculate axes
}

result_df$x <- x.v
result_df$y <- y.v 


order.v <- rownames(merged.seuratobj.sct@meta.data)
colnames(UMAP_RNA) <- c("x","y")   
UMAP_RNA$x <- NA
UMAP_RNA$y <- NA

# MERFISH spatial loci substitution
UMAP_RNA$x[which(is.na(str_match(colnames(dist_matrix_mat),"TSO")))] <- loc.df[rownames(UMAP_RNA)[which(is.na(str_match(colnames(dist_matrix_mat),"TSO")))],'x']
UMAP_RNA$y[which(is.na(str_match(colnames(dist_matrix_mat),"TSO")))] <- loc.df[rownames(UMAP_RNA)[which(is.na(str_match(colnames(dist_matrix_mat),"TSO")))],'y']

# our RNA spatial loci substitution
rownames(result_df) <- result_df$ide
UMAP_RNA$x[which(!is.na(str_match(rownames(UMAP_RNA),"TSO")))] <- result_df[rownames(UMAP_RNA)[which(!is.na(str_match(rownames(UMAP_RNA),"TSO")))],'x']
UMAP_RNA$y[which(!is.na(str_match(rownames(UMAP_RNA),"TSO")))] <- result_df[rownames(UMAP_RNA)[which(!is.na(str_match(rownames(UMAP_RNA),"TSO")))],'y']


UMAP_RNA <- UMAP_RNA[order.v,]

colnames(UMAP_RNA) <- c("umap_1","umap_2")
UMAP_RNA <- as.matrix(UMAP_RNA)
# replace to umap dimension
merged.seuratobj.sct@reductions$umap@cell.embeddings <- UMAP_RNA 

saveRDS(merged.seuratobj.sct,file=paste0(args$neuron_type,"_obj.label_transfer_",args$label,".loci_transfer.the_nearst_1_cell.rds"))

