#### 01-job.sh: 
This shell script is a SLURM batch submission script designed to run the 01-scRNAseq.R pipeline.

#### 01-scRNAseq.R
This script processes single-cell RNA-seq data from the Joint Cabernet project by reading filtered count matrices and metadata, then creating Seurat objects for downstream analysis. It supports batch processing via a metadata file or single-sample analysis through command-line input. The pipeline filters cells based on provided RNA statistics, matches gene identifiers to a mouse reference, and generates quality control plots including feature counts and mitochondrial content.

#### 02-integration.R
This script performs integration of multiple single-cell RNA-seq Seurat objects using the SCT normalization and reciprocal PCA (RPCA) approach. It preprocesses each dataset with SCTransform, identifies integration anchors, and generates an integrated Seurat object. The script then conducts dimensionality reduction (PCA, t-SNE, UMAP), clustering, and marker gene identification for each cluster.

#### 02-integration.sh
This bash script is a SLURM job submission script designed to run the 02-integration.R analysis pipeline.