#### 01-job.sh: 
This shell script is a SLURM batch submission script designed to run the 01-scRNAseq.R pipeline.

#### 01-scRNAseq.R
This script processes single-cell RNA-seq data from the Joint Cabernet project by reading raw count matrices, then filtering with RNA QC, and finally creating Seurat objects and filtering by nfeature and ncounts for downstream analysis. It supports batch processing via a metadata file or single-sample analysis through command-line input. The pipeline filters cells based on provided RNA statistics, matches gene identifiers to a mouse reference, and generates quality control plots including feature counts and mitochondrial content.

#### 01-meta.txt
input of 01-job.sh

#### 02-change_QC_before_integration.R
The code modified the quality control for RNA by removing cells that did not meet the nfeature and ncount filtering criteria.

#### 02-change_QC_before_integration.sh
This bash script is a SLURM job submission script designed to run the 02-change_QC_before_integration.R analysis pipeline.