01-preprocessing.r  :  Preprocessing the data. Generating our young and aged RNA data after integration with zeng and filter out low-quality genes.
The filter condition of genes is that the count value is greater than 500 and the number of cells with cpm value greater than 0.5 is greater than 2.

02-dividing_metacells.r  :  Dividing single cells into metacells.
Note：If an error occurs when loading the hdWGCNA Package, for example, Package 'Rcpp' version 1.0.11 cannot be unloaded, re-open an R window.

03-screen_aged_DEG.r  :  Screening aging related differential genes by DESeq2.