## options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)
print(args)
metafile <- args[1]
species <- args[2]
stat <- args[3]

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

EnsDb <- c(human="GRCh38", mouse="mm10", monkey="macFas5", wheat='wheat', arabidopsis="tair10")
commonSpecies <- names(EnsDb)
if(species %in% commonSpecies){
  cat("Species ", species, "was selected for analysis\n")
}else{
  cat("Currently supported species: ", paste(commonSpecies, sep="", collapse = ","), "\n")
  q()
}

dir.create("rds", showWarnings = F)
dir.create("plots", showWarnings = F)

runAnalysis <- function(meta){
  if(length(meta) == 1){
    datasetID <- datasetName <- meta[1]
  }else{
    datasetName <- meta[1]
    datasetID <- meta[2]
  }
  if(file.exists(sprintf("rds/%s_seurat.rds", datasetID))){
    print(sprintf("Result file rds/%s_seurat.rds exists. Skipping ...", datasetID))
    return(NULL)
  }
  
  cat("Running analysis for dataset ", datasetName, "\n")
  
  ## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(Seurat)
  
  mat <- readRDS(sprintf("input/%s.rds",datasetName))
  mat <- as(as.matrix(mat), "dgCMatrix")

  QC_stat <- read.csv(stat,header=T)
  QC_stat_filter <- QC_stat %>% dplyr::filter(QC == 1)
  mat <- mat[,colnames(mat) %in% QC_stat_filter$SampleID]
 
  mm10_geneID = read.delim("../../../03.data/01.ref/mm10.genes.bed",header=F) 
  colnames(mm10_geneID) = c("chr","start","end","gene_id","gene_name","gene_type")
  
  
  ##### Create object #####
  seuratObj <- CreateSeuratObject(counts=mat, project=datasetName, min.cells=0, min.features=0)
  seuratObj <- RenameCells(object = seuratObj, add.cell.id = paste0(datasetID,"@@"))
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(
    seuratObj, 
    features = mm10_geneID[mm10_geneID$chr == "chrM","gene_id"],
    assay = "RNA"
  )
  pdf(file=sprintf("plots/%s_qc_plots.pdf", datasetID), width=6, height=8.27)
  p1 <- VlnPlot(seuratObj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  p2 <- FeatureScatter(seuratObj, feature1="nCount_RNA", feature2="nFeature_RNA")
  print(p1)
  print(p2)
  dev.off()

  ### filter 
  # 3. calculate IQR threshold
  calculate_thresholds <- function(metric) {
      q <- quantile(seuratObj[[metric]][, 1], probs = c(0.25, 0.75))
      iqr <- IQR(seuratObj[[metric]][, 1])
      lower <- q[1] - 1.5 * iqr
      upper <- q[2] + 1.5 * iqr
      return(c(lower, upper))
  }

  nFeature_thresholds <- calculate_thresholds("nFeature_RNA")
  nCount_thresholds <- calculate_thresholds("nCount_RNA")

  # 4. filter cells
  seuratObj <- subset(seuratObj,
      subset = nFeature_RNA > nFeature_thresholds[1] & 
           nFeature_RNA < nFeature_thresholds[2] &
           nCount_RNA > nCount_thresholds[1] & 
           nCount_RNA < nCount_thresholds[2] 
  )
  print(paste0("the number of retain cells is ",dim(seuratObj@meta.data)[1]))
  retain_cell = colnames(seuratObj@assays$RNA$counts)
  retain_gene = rownames(seuratObj@assays$RNA$counts)
  mm10_geneID = mm10_geneID[-which(mm10_geneID$chr %in% c("chrX","chrY","chrM")),]
  mm10_geneID$geneid_without_version = unlist(lapply(mm10_geneID$gene_id, function(x) strsplit(x,'\\.')[[1]][1]))
  feature.names <- mm10_geneID[,c(3:4)]
  refGenes <- feature.names[,2]
  names(refGenes) <- feature.names[,1]
  refGenes[refGenes==""] <- names(refGenes)[refGenes==""]
  refGenes[duplicated(refGenes)] <- names(refGenes[duplicated(refGenes)])
  saveRDS(refGenes, sprintf("rds/%s_genes.rds", datasetID))

  mat = mat[intersect(mm10_geneID$gene_id,retain_gene),]
  mat = mat[,unlist(lapply(retain_cell, function(x) strsplit(x,'@@_')[[1]][2]))]
  toppc <- 30
  seuratObj <- CreateSeuratObject(counts=mat, project=datasetName, min.cells=0, min.features=0)
  seuratObj <- RenameCells(object = seuratObj, add.cell.id = paste0(datasetID,"@@"))
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt-")
  seuratObj <- SCTransform(seuratObj, vars.to.regress="percent.mt", verbose = FALSE)
  seuratObj <- RunPCA(seuratObj, verbose = FALSE)
  seuratObj <- RunUMAP(seuratObj, dims=1:toppc, n.components=3, verbose = FALSE)
  seuratObj <- RunTSNE(seuratObj, dims=1:toppc, dim.embed=3, verbose=FALSE)
  seuratObj <- FindNeighbors(seuratObj, dims = 1:toppc, verbose = FALSE)
  seuratObj <- FindClusters(seuratObj, verbose = FALSE) 
  
  saveRDS(seuratObj, file=sprintf("rds/%s_seurat.rds", datasetID))
}

## mode 1
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
}else{ ## mode 2
  paras <- unlist(strsplit(metafile, split=":"))
  if(length(paras)<2){
    print("Wrong format of metadata [e.g., 'sample:file']")
    q()
  }
  runAnalysis(paras[1:2])
}
