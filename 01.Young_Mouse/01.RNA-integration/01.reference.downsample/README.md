#### 01-generate_RNA_expression_counts_from_h5ad.R: 
This script reads multiple .h5ad files using Python's Scanpy through the reticulate package in R, extracts and transposes RNA expression count matrices, and saves them as .rds files. It constructs unique cell names by combining barcodes with library labels and is designed for efficient integration with downstream R-based single-cell analysis workflows.

#### 02-generate_metaData_from_h5ad.R:
This script uses the reticulate package to interface with Python's Scanpy and extract cell-level metadata (.obs) from multiple .h5ad files. For each dataset, it creates a unique cell identifier by concatenating the cell barcode and library label, then writes the full metadata table to a tab-delimited .txt file. 

#### 03-generate_downsampled_metadata.ipynb:
This notebook processes and integrates metadata from multiple single-cell RNA-seq .h5ad files to generate a curated and downsampled metadata table. It merges region-specific cell metadata with global annotations (including anatomical labels, region labels, class/subclass identities), filters for specific brain regions of interest, and selects subclasses with sufficient cell counts (n > 100). A stratified random sampling strategy is then applied to select up to 1000 cells per subclass-region combination. The final output is saved as an .rds file for downstream analysis.

#### 04-generate_downsampled_RNA_expression_counts.ipynb:
This notebook subsets full RNA expression matrices to include only the downsampled cells selected in the curated metadata file. For each brain region, it loads the original expression matrix, filters the columns to retain only selected cells, and saves the resulting matrix as an .rds file.
