#### 01-integration.R: 
This script performs integrated analysis of multiple single-cell RNA-seq datasets using Seurat’s SCT normalization and RPCA-based integration workflow. It reads preprocessed Seurat objects, identifies shared variable features, computes integration anchors, and merges datasets into a single integrated object. The pipeline then runs dimensionality reduction (PCA, t-SNE, UMAP), clustering, and identifies cluster-specific marker genes, saving intermediate and final Seurat objects along with marker gene lists for downstream interpretation.

#### 01-integration.sh
This shell script is a SLURM job submission script designed to run an R integration analysis (01-integration.R) with specified metadata (01-meta.txt).

#### 01-meta.txt
input of 01-integration.sh

#### 02-label_transfer.R
The code assigns the major class and subclass labels from the reference data with the largest proportion within the same Seurat cluster to the Joint Cabernet cells in that cluster. We refer to this step as label transfer. Its main purpose is to annotate the Joint Cabernet cells.

#### 02-label_transfer.sh
This shell script is a SLURM job submission script designed to run 02-label_transfer.R