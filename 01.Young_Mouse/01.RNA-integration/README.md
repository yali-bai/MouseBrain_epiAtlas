#### 01.reference.downsample
This folder contains scripts and notebooks for processing single-cell RNA-seq data starting from raw .h5ad files. It extracts RNA expression count matrices and cell-level metadata via Scanpy in R, generates unique cell identifiers, integrates and downsamples the data, and outputs curated expression matrices and metadata files suitable for downstream single-cell analysis in R.

#### 02.Zeng_10x.RNA
This folder includes an R script for quality control and preprocessing of single-cell RNA-seq data using the Seurat framework. It supports both batch and single-sample processing, filters gene IDs based on mouse reference, standardizes gene annotations, creates Seurat objects, and produces basic QC plots to facilitate further analysis.

#### 03.Joint-Cabernet.RNA
This folder contains analysis scripts and SLURM job scripts for processing Joint Cabernet single-cell RNA-seq data. The pipeline includes data loading, filtering, Seurat object creation, QC visualization, and multi-dataset integration using SCT normalization and RPCA. It performs dimensionality reduction, clustering, and marker gene identification.

#### 04.Joint-Cabernet.Zeng_10X_RNA.integration
This folder focuses on integrated analysis of multiple single-cell RNA-seq datasets using Seurat’s SCT normalization and RPCA-based integration workflow. It includes R scripts and SLURM submission scripts for merging datasets, performing dimensionality reduction, clustering, and marker gene detection.