# "indir" is a custom input path, and "outdir" is a custom output path.
indir="" # h5ad file path
outdir="./"

##########   1. compute_5mC
#We have 5hmC and 5hmC+5mC sequencing data, and we compute 5mC data by '5hmC+5mC' -5hmC
import scanpy as sc
import pandas as pd  
  
def compute_5mC(var_dim,mc_type):  
    #Load two h5ad files
    hmC = sc.read_h5ad(f'{indir}/5hmC_{var_dim}.{mc_type}.pass_total_QC.h5ad')  
    mC = sc.read_h5ad(f'{indir}/5mC_{var_dim}.{mc_type}.pass_total_QC.h5ad')  
    #Check whether row names and column names are consistent
    print(hmC.obs_names.equals(mC.obs_names))
    print(hmC.var_names.equals(mC.var_names))
    info=pd.read_csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv")
    info.index='allc_'+info['hmC_SampleID']
    hmC.obs_names=info['unique_id'][hmC.obs_names]
    info.index='allc_'+info['mC_SampleID']
    mC.obs_names=info['unique_id'][mC.obs_names]    
    mC=mC[hmC.obs_names,:]
    print(hmC.obs_names.equals(mC.obs_names))
    #compute true5mC
    true5mC_matrix = mC.X - hmC.X
    true5mC_matrix[true5mC_matrix<0]=0
    true5mC_adata = sc.AnnData(X=true5mC_matrix) 
    true5mC_adata.obs.index=hmC.obs.index
    true5mC_adata.var.index=hmC.var.index
    #add obs and var
    info = info[info['total_QC'] == 1]
    info.index=info['unique_id']
    info = info.reindex(hmC.obs_names)

    columns_to_add = ['RNA_SampleID','hmC_SampleID','mC_SampleID','unique_id','Region','Brain_Region','class_label','subclass_label','Neuron_non_neuron','three_class']  
    selected_columns = info[columns_to_add]
    true5mC_adata.obs[columns_to_add] = selected_columns  
    true5mC_adata.obs['Library']='true-5mC'
    true5mC_adata.var = hmC.var
    true5mC_adata.write_h5ad(f'{outdir}/true_5mC_{var_dim}.{mc_type}.pass_total_QC.h5ad')

    hmC.obs[columns_to_add] = selected_columns
    hmC.obs['Library']='5hmC'
    mC.obs[columns_to_add] = selected_columns
    mC.obs['Library']='5mC'
    hmC.write_h5ad(f'{outdir}/5hmC_{var_dim}.{mc_type}.pass_total_QC.h5ad')
    mC.write_h5ad(f'{outdir}/5mC_{var_dim}.{mc_type}.pass_total_QC.h5ad')


var_dim='genebody'
mc_types = ['CH'] #,'CH'
# Loops through all parameter combinations and calls the function
for mc_type in mc_types:
    print(f"Running for mc_type: {mc_type}, var_dim: {var_dim}")
    compute_5mC(
        mc_type=mc_type,
        var_dim=var_dim                
    )