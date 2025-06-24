# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

################################  shuffle correlation   #################
#Randomly perturbed DNA sample names 5000 times to calculate the correlation between gene RNA expression and DNA methylation in the random state.


###################    1.random sort DNA samples   
## shuffled_correlation is prepared here to save the sequence after 5000 random perturbations of DNA samples for subsequent extraction

import scanpy as sc
import pandas as pd
from scipy.sparse import issparse
import re 
import pickle


def random_sort(datatype,var_dim,mc_type):
    #RNA
    Joint_Cabernet_RNA=sc.read_h5ad(f'{indir}/Joint_Cabernet_RNA_annotated_latest.h5ad')
    if issparse(Joint_Cabernet_RNA.X):  
        X_dense = Joint_Cabernet_RNA.X.toarray()  
    else:  
        X_dense = Joint_Cabernet_RNA.X
    Joint_Cabernet_RNA_count= pd.DataFrame(X_dense, columns=Joint_Cabernet_RNA.var_names) 
    Joint_Cabernet_RNA_count.index=Joint_Cabernet_RNA.obs.index
    #Extract neuron cells
    Joint_Cabernet_RNA_count=Joint_Cabernet_RNA_count.loc[Joint_Cabernet_RNA.obs.index[Joint_Cabernet_RNA.obs['Neuron_non_neuron']=="Neuron"],:]
    Joint_Cabernet_RNA_count.index=[re.sub(r'.*@@', 'allc', s) for s in Joint_Cabernet_RNA_count.index]
    #Extract paired QC
    pairedQC=pd.read_csv('../../../04.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv')
    pairedQC.index='allc_' + pairedQC['RNA'].astype(str) 
    common=pairedQC[pairedQC['total_QC'] == 1].index.intersection(Joint_Cabernet_RNA_count.index)
    RNA=Joint_Cabernet_RNA_count.loc[common,:]
    RNA.index=pairedQC['Unique_ID_match'][common]
    #Filter gene
    col_sums = RNA.sum()  
    RNA = RNA.loc[:, col_sums != 0] # Filter out columns whose sum is not 0
    zero_proportions = RNA.apply(lambda col: (col == 0).sum() / len(col))  #Calculate the proportion of 0
    RNA = RNA.loc[:, zero_proportions < 0.8]

    #mC
    mc=sc.read_h5ad(f'{indir}/{datatype}_{var_dim}.{mc_type}.pass_paired_QC.h5ad')
    if issparse(mc.X):  
        X_dense = mc.X.toarray()  
    else:  
        X_dense = mc.X
    DNA= pd.DataFrame(X_dense, columns=mc.var_names)
    DNA.index=mc.obs['Unique_ID']
    zero_proportions2 = DNA.isna().sum() / len(DNA)  #Calculate the proportion of NA
    DNA = DNA.loc[:, zero_proportions2 < 0.8]
    #Extracts the common row name column name
    common_rows = RNA.index.intersection(DNA.index) 
    common_columns = RNA.columns.intersection(DNA.columns)
    # Sort and filter data boxes by common row and column names
    RNA_common = RNA.loc[common_rows, common_columns]
    DNA_common = DNA.loc[common_rows, common_columns]
    DNA_common['OriginalIndex']=DNA_common.index
    # Initializes a list to store the index order after each sort
    sorted_indices = [] 
    # Do 5000 random sorting  
    for _ in range(5000):  
        # Randomly sort the DataFrame, but keep the 'OriginalIndex' column for tracking
        shuffled_df = DNA_common.sample(frac=1).reset_index(drop=True)  
        
        # Record the sorted index order (here we are actually recording the order of the 'OriginalIndex' column)  
        sorted_order = shuffled_df['OriginalIndex'].tolist()  
        sorted_indices.append(sorted_order)

    # pickle data to a file  
    path=f"{outdir}/random_sort"
    with open(f'{path}/{datatype}_{var_dim}_{mc_type}.pkl', 'wb') as file:  
        pickle.dump(sorted_indices, file)  
    # Note: pickle files are binary and are not recommended for data exchange as it depends on Python version and platform

datatypes=['5hmC','5mC','true_5mC']
var_dims = ['genebody']
mc_types = ['CG','CH']
for datatype in datatypes:  
    for mc_type in mc_types:
        for var_dim in var_dims:
            print(f"Running for datatype: {datatype}, mc_type: {mc_type}, var_dim: {var_dim}")
            random_sort(
                datatype=datatype,
                var_dim=var_dim, 
                mc_type=mc_type               
            )




#########################   2. shuffled correlation
# Here 5000 random perturbations are divided into 100 times to execute, each time extracting 50 out of 5000 sample sequences
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse
import re 
from scipy.stats import pearsonr
import pickle


def allcell_shuffle_correlation(datatype,var_dim,mc_type,num,n1,n2):
    #RNA
    Joint_Cabernet_RNA=sc.read_h5ad(f'{indir}/Joint_Cabernet_RNA_annotated_latest.h5ad')
    if issparse(Joint_Cabernet_RNA.X):  
        X_dense = Joint_Cabernet_RNA.X.toarray()  
    else:  
        X_dense = Joint_Cabernet_RNA.X
    Joint_Cabernet_RNA_count= pd.DataFrame(X_dense, columns=Joint_Cabernet_RNA.var_names) 
    Joint_Cabernet_RNA_count.index=Joint_Cabernet_RNA.obs.index
    #Extract neuron cells
    Joint_Cabernet_RNA_count=Joint_Cabernet_RNA_count.loc[Joint_Cabernet_RNA.obs.index[Joint_Cabernet_RNA.obs['Neuron_non_neuron']=="Neuron"],:]
    Joint_Cabernet_RNA_count.index=[re.sub(r'.*@@', 'allc', s) for s in Joint_Cabernet_RNA_count.index]
    #Extract paired QC
    pairedQC=pd.read_csv('../../../04.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv')
    pairedQC.index='allc_' + pairedQC['RNA'].astype(str) 
    common=pairedQC[pairedQC['total_QC'] == 1].index.intersection(Joint_Cabernet_RNA_count.index)
    RNA=Joint_Cabernet_RNA_count.loc[common,:]
    RNA.index=pairedQC['Unique_ID_match'][common]
    #Filter gene
    col_sums = RNA.sum()  
    RNA = RNA.loc[:, col_sums != 0] # Filter out columns whose sum is not 0 
    zero_proportions = RNA.apply(lambda col: (col == 0).sum() / len(col))  #Calculate the proportion of 0
    RNA = RNA.loc[:, zero_proportions < 0.8]

    #DNA 
    mc=sc.read_h5ad(f'{indir}/{datatype}_{var_dim}.{mc_type}.pass_paired_QC.h5ad')
    if issparse(mc.X):  
        X_dense = mc.X.toarray()  
    else:  
        X_dense = mc.X
    DNA= pd.DataFrame(X_dense, columns=mc.var_names)
    DNA.index=mc.obs['Unique_ID']
    zero_proportions2 = DNA.isna().sum() / len(DNA)  #Calculate the proportion of NA
    DNA = DNA.loc[:, zero_proportions2 < 0.8]
    #Extracts the common row name column name
    common_rows = RNA.index.intersection(DNA.index) 
    common_columns = RNA.columns.intersection(DNA.columns)
    # Sort and filter data boxes by common row and column names
    RNA_common = RNA.loc[common_rows, common_columns]
    DNA_common = DNA.loc[common_rows, common_columns]


    ###############Random disturbance
    path=f"{outdir}/random_sort"
    with open(f'{path}/{datatype}_{var_dim}_{mc_type}.pkl', 'rb') as file:  
        sorted_indices = pickle.load(file)
    # Initializes a list to store the correlation results of all iterations
    all_correlations = []
    # Enter each start and end
    n1=int(n1)
    n2=int(n2)
    for i in range(n1,n2 + 1): 
        print(i) 
        if i==5000:
            i=0
        shuffled_dna = DNA_common.loc[sorted_indices[i],:]      # Scramble the rows of the DNA matrix
        # Initializes a dictionary to store the dependencies of the current iteration 
        current_correlations = {}  
        # Go through all genes 
        for gene in RNA_common.columns:  
            # Extract expression data of corresponding genes in RNA and DNA
            rna_vals = RNA_common[gene]  
            dna_vals = shuffled_dna[gene] 
            dna_vals.index=rna_vals.index
            df2=pd.DataFrame({
            'DNA':dna_vals,
            'RNA':rna_vals
            })
            df_cleaned2 = df2.dropna(subset=['DNA'])
            # Calculate Pearson correlation coefficient 
            corr, _ = pearsonr(df_cleaned2['DNA'], df_cleaned2['RNA'])  
            # Store results  
            current_correlations[gene] = corr  
        # Adds the correlation results of the current iteration to the list  
        all_correlations.append(current_correlations)  
  
    gene_ids = list(all_correlations[0].keys())
    dfs = [pd.DataFrame({gene: [d[gene] for d in all_correlations]}) for gene in gene_ids] 
    # Merge these dataframes using concat, axis=1 means merge along the column direction
    result_df = pd.concat(dfs, axis=1) 
    result_df=result_df.T
    result_df.columns=[f'cor{j}' for j in range(n1, n2 + 1)]
    # The results are stored 100 times each and then consolidated
    result_df.to_csv(f'{outdir}/shuffled_{datatype}_{mc_type}_{var_dim}_gene_correlation_results({num}).csv')


datatypes=['5hmC','5mC','true_5mC']
var_dims = ['genebody']
mc_types = ['CG','CH']
total = 5000
parts = 100
per_part = total // parts  
for datatype in datatypes:  
    for mc_type in mc_types:
        for var_dim in var_dims:
            for i in range(parts):  # Generates a sequence of integers from 0 to parts-1
                #Python's index starts at 0, but here we are calculating the actual value
                first = i * per_part + 1
                # If it is the last copy, you need to deal with the remaining balance that may exist
                if i == parts - 1:
                    # If it is the last copy, the last number is the total
                    last = total
                else:
                    # Otherwise, the last number is the last position of the current number
                    last = first + per_part - 1
                num = i + 1
                print(f"Part {num}: First = {first}, Last = {last}")
                print(f"Running for datatype: {datatype}, mc_type: {mc_type}, var_dim: {var_dim}")
                allcell_shuffle_correlation(datatype,var_dim,mc_type,num,first,last)

# The extraction of each sample name starts from 1 when the loop is performed here, but python is 0-based,
# Therefore, i=0 when i=5000 is added to the allcell_shuffle_correlation function, and the missing i=0 is replaced








