#########    All "our" in the following code refers to Joint Cabernet.
########RNA expression and DNA methylation dotplot plot data preparation in each subclass of gene AC132685.1.
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse
import re 

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

#gene_list=['AC132685.1']
gene='ENSMUSG00000068151.7'
#RNA
our_RNA=sc.read_h5ad(f'{indir}/our_RNA_annotated_latest.h5ad')
if issparse(our_RNA.X):  
    X_dense = our_RNA.X.toarray()  
else:  
    X_dense = our_RNA.X
our_RNA_count= pd.DataFrame(X_dense, columns=our_RNA.var_names) 
our_RNA_count.index=our_RNA.obs.index
# Extract neuron cells
our_RNA_count=our_RNA_count.loc[our_RNA.obs.index[our_RNA.obs['Neuron_non_neuron']=="Neuron"],:]
our_RNA_count.index=[re.sub(r'.*@@', 'allc', s) for s in our_RNA_count.index]
# paired QC was extracted
pairedQC=pd.read_csv('../../../input/01-youth/RNA_DNA_match_name_QC_class_label.csv')
pairedQC.index='allc_' + pairedQC['RNA'].astype(str) 
common=pairedQC[pairedQC['total_QC'] == 1].index.intersection(our_RNA_count.index)
RNA=our_RNA_count.loc[common,:]
RNA.index=pairedQC['Unique_ID_match'][common]
col_sums = RNA.sum()  
#Filter out columns whose sum is not 0
RNA = RNA.loc[:, col_sums != 0] 


datatypes=['5hmC','5mC','true_5mC']
var_dim="genebody"
mc_type="CG"
result_df=pd.DataFrame()
for datatype in datatypes:  
    #DNA 
    mc=sc.read_h5ad(f'{indir}/{datatype}_{var_dim}.{mc_type}.pass_paired_QC.h5ad')
    if issparse(mc.X):  
        X_dense = mc.X.toarray()  
    else:  
        X_dense = mc.X
    DNA= pd.DataFrame(X_dense, columns=mc.var_names)
    DNA.index=mc.obs['Unique_ID']
    #Find the common row and column names
    common_rows = RNA.index.intersection(DNA.index) 
    common_columns = RNA.columns.intersection(DNA.columns)
    # Sort and filter data boxes by common row and column names
    RNA_common = RNA.loc[common_rows, common_columns]
    DNA_common = DNA.loc[common_rows, common_columns]
    pairedQC_1=pairedQC[pairedQC['total_QC'] == 1]
    pairedQC_1.index=pairedQC_1['Unique_ID_match']
    pairedQC_1=pairedQC_1.loc[common_rows,:]
    df=pd.DataFrame({
        'gene_id':gene,
        'DNA':DNA_common.loc[:,gene],
        'RNA':RNA_common.loc[:,gene],
        'datatype':datatype,
        'subclass':pairedQC_1['subclass_label'],
        'class':pairedQC_1['class_label'],
        'three_class':pairedQC_1['three_class_label']
        })
    df_cleaned = df.dropna(subset=['DNA'])
    result_df=result_df.append(df_cleaned)
    result_df.to_csv('../../../output/01-youth/02-correlation_calculation/genebody_CG_ENSMUSG00000068151.7.csv')
    print('genebody_CG_ENSMUSG00000068151.7:finished')


