#### 01-integration.R: 
This script performs integrated analysis of multiple single-cell RNA-seq datasets using Seurat’s SCT normalization and RPCA-based integration workflow. It reads preprocessed Seurat objects, identifies shared variable features, computes integration anchors, and merges datasets into a single integrated object. The pipeline then runs dimensionality reduction (PCA, t-SNE, UMAP), clustering, and identifies cluster-specific marker genes, saving intermediate and final Seurat objects along with marker gene lists for downstream interpretation.

#### 01-integration.sh
This shell script is a SLURM job submission script designed to run an R integration analysis (01-integration.R) with specified metadata (01-meta.txt).