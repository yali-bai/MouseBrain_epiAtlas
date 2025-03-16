#########    All "our" in the following code refers to Joint Cabernet.

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

## options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)
print(args)
metafile <- args[1]
species <- args[2]

usage <- function(){
  cat("\n***************************************************************\n")
  cat("* eg 1: scRNAseq.R [OPTIONS]
                            metafile species[human/mouse]\n")
  cat("* eg 2: scRNAseq.R [OPTIONS]
                            'sample:file' species[human/mouse]\n")
  cat("***************************************************************\n")
}

if(length(args)<2){
  print("[ERROR] Less than two arguments ...", quote = F)
  usage()
  q()
}

dir.create(paste0(outdir,"/rds/"), showWarnings = F)
dir.create(paste0(outdir,"/plots/"), showWarnings = F)

runAnalysis <- function(meta){
  if(length(meta) == 1){
    datasetID <- datasetName <- meta[1]
  }else{
    datasetName <- meta[1]
    datasetID <- meta[2]
  }
  if(file.exists(sprintf("%s/rds/%s_seurat.rds",outdir, datasetID))){
    print(sprintf("Result file %s/rds/%s_seurat.rds exists. Skipping ...",outdir, datasetID))
    return(NULL)
  }
  
  cat("Running analysis for dataset ", datasetName, "\n")
  
  ## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(data.table)
  
  mat <- readRDS(sprintf("%s/%s_counts.rds",indir,datasetName)) 
  mat <- as(as.matrix(mat), "dgCMatrix")
  #mouse_geneID <- read.delim("../data/mm10.genes_duplicated.bed",header=T)
  #mat <- mat[rownames(mat) %in% mouse_geneID$gene_name,]
  #mat <- mat[!duplicated(rownames(mat)),]
  #rownames(mat) <- mouse_geneID$gene_id[match(rownames(mat),mouse_geneID$gene_name)]
  RNA_stat <- read.delim(sprintf("%s/RNA_stat_filter.txt",indir))  
  mat <- mat[,colnames(mat) %in% RNA_stat$SampleID]
  mm10_geneID <- read.delim("../../../input/reference_genome/mm10_vM18.genes.bed",header=T)  
  feature.names <- mm10_geneID
  refGenes <- feature.names$gene_name
  names(refGenes) <- feature.names$gene_id
  refGenes[refGenes==""] <- names(refGenes)[refGenes==""]
  refGenes[duplicated(refGenes)] <- names(refGenes[duplicated(refGenes)])
  saveRDS(refGenes, sprintf("%s/rds/%s_genes.rds",outdir, datasetID))  
  
  ##### Create object #####
  library(Seurat)
  seuratObj <- CreateSeuratObject(counts=mat, project=datasetName, min.cells=0, min.features=0)
  seuratObj <- RenameCells(object = seuratObj, add.cell.id = paste0(datasetID,"@@"))
  pdf(file=sprintf("%s/plots/%s_qc_plots.pdf",outdir, datasetID), width=6, height=8.27)  
  p1 <- VlnPlot(seuratObj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  p2 <- FeatureScatter(seuratObj, feature1="nCount_RNA", feature2="nFeature_RNA")
  print(p1)
  print(p2)
  dev.off()
  
  saveRDS(seuratObj, file=sprintf("%s/rds/%s_seurat.rds",outdir, datasetID))  
}

if(file.exists(metafile)){
  metadata <- read.delim(metafile, header = F)
  rownames(metadata) <- metadata[,1]
  if(ncol(metadata)<2){
    print("Wrong format of metafile")
    q()
  }
  library(doMC)
  registerDoMC(30)
  mclapply(1:nrow(metadata), function(i){
    runAnalysis(unlist(metadata[i,1:2]))
  }, mc.cores = 10) 
}else{
  paras <- unlist(strsplit(metafile, split=":"))
  if(length(paras)<2){
    print("Wrong format of metadata [e.g., 'sample:file']")
    q()
  }
  runAnalysis(paras[1:2])
}



