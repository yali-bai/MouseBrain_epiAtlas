# "indir" is a custom input path, and "outdir" is a custom output path.
indir="./"
# outdir=""

################    2.  DNA-RNA correlation of genes in all cells 
# Pearson correlation coefficients and P-values of RNA expression and DNA methylation of each gene in all cells with 5hmC/5hmC+5mC/5mC, genebody, CG/CH data were calculated.
import sys
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse
import re 
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

import anndata as ad
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import numpy as np
import os

import anndata as ad
import pandas as pd
import scanpy as sc
import scglue

pandas2ri.activate()


def allcell_correlation(datatype,var_dim,mc_type):
    #RNA
    rds_data = ro.r['readRDS']("sc_RNA_expression_log1p_normalized_data.rds") 
    ro.globalenv['dataframe'] = rds_data
    rds_data_temp= pandas2ri.rpy2py(rds_data)
    Joint_Cabernet_RNA=ad.AnnData(rds_data_temp.T)

    if issparse(Joint_Cabernet_RNA.X):  
        X_dense = Joint_Cabernet_RNA.X.toarray()  
    else:  
        X_dense = Joint_Cabernet_RNA.X
    Joint_Cabernet_RNA_count= pd.DataFrame(X_dense, columns=Joint_Cabernet_RNA.var_names) 
    Joint_Cabernet_RNA_count.index=Joint_Cabernet_RNA.obs.index

    pairedQC=pd.read_csv('../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv')
    pairedQC.index=pairedQC['RNA_SampleID'].astype(str) 
    common=pairedQC[pairedQC['total_QC'] == 1].index.intersection(Joint_Cabernet_RNA_count.index)
    RNA=Joint_Cabernet_RNA_count.loc[common,:]
    RNA.index=pairedQC['unique_id'][common]
    #Filter gene
    col_sums = RNA.sum()  
    RNA = RNA.loc[:, col_sums != 0] # Filter out columns whose sum is not 0
    zero_proportions = RNA.apply(lambda col: (col == 0).sum() / len(col))  #Calculate the proportion of 0
    RNA = RNA.loc[:, zero_proportions < 0.8]

    #DNA 
    mc=sc.read_h5ad(f'{indir}/{datatype}_{var_dim}.{mc_type}.pass_total_QC.h5ad')
    if issparse(mc.X):  
        X_dense = mc.X.toarray()  
    else:  
        X_dense = mc.X
    DNA= pd.DataFrame(X_dense, columns=mc.var_names)
    DNA.index=mc.obs['unique_id']
    zero_proportions2 = DNA.isna().sum() / len(DNA)  #Calculate the proportion of NA
    DNA = DNA.loc[:, zero_proportions2 < 0.8]
    #Extracts the common row name column name
    common_rows = RNA.index.intersection(DNA.index) 
    common_columns = RNA.columns.intersection(DNA.columns)
    # Sort and filter data boxes by common row and column names
    RNA_common = RNA.loc[common_rows, common_columns]
    DNA_common = DNA.loc[common_rows, common_columns]
    #Create a data frame to store the results
    correlation_results_all = pd.DataFrame(index=RNA_common.columns, columns=["datatype","var_dim","mc_type","paired_cell_num",'Correlation', 'P-value'])
    # Calculate the correlation of each gene
    for gene in RNA_common.columns:
        df=pd.DataFrame({
            'DNA':DNA_common.loc[:,gene],
            'RNA':RNA_common.loc[:,gene]
            })
        df_cleaned = df.dropna(subset=['DNA'])
        #all
        if len(df_cleaned.index) > 2:  # Ensure that there are enough data points for correlation calculations
            paired_cell_num=len(df_cleaned.index)
            correlation, p_value = pearsonr(df_cleaned['DNA'], df_cleaned['RNA'])
        else:
            paired_cell_num,correlation, p_value = len(df_cleaned.index), float('nan'), float('nan')  
        correlation_results_all.loc[gene] = [datatype,var_dim,mc_type,paired_cell_num,correlation, p_value]
        print(f'{gene} finished')
    # Gets all the original p-values
    p_values = correlation_results_all['P-value'].values
    # Multiple comparison correction was performed using the Benjamini-Hochberg method
    _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    # Adds the corrected p-value to the result DataFrame
    correlation_results_all['Adjusted P-value'] = corrected_p_values
    # Save the results to a CSV file
    correlation_results_all.to_csv(f'all_cell_correlation/all_cells_{datatype}_{mc_type}_{var_dim}_gene_correlation_results.csv')



#datatypes=['5hmC','5mC','true_5mC']
#var_dims = ['genebody']
#mc_types = ['CG','CH']
if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)

    datatype = args[0]  # "5hmC"
    var_dim = args[1] # "geneslop2k"
    mc_type = args[2] # "CG"

    

    allcell_correlation(datatype,var_dim,mc_type)
