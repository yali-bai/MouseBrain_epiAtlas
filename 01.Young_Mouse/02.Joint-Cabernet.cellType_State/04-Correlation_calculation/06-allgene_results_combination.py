##########################    all gene results combination
import pandas as pd 

datatypes=['5hmC','5mC','true_5mC']
var_dims = ['genebody']
mc_types = ['CG','CH']
df=pd.DataFrame()
for mc_type in mc_types:
    for datatype in datatypes:
        for var_dim in var_dims:
            print(f"Running for mc_type: {mc_type}, datatype: {datatype}, var_dim: {var_dim}")
            data=pd.read_csv(f"../../../output/01-youth/02-correlation_calculation/shuffled_correlation/shuffled_{datatype}_{mc_type}_{var_dim}_gene_correlation_total_results.csv")
            df_filtered = data.dropna(subset=['Correlation'])
            df=df.append(df_filtered) 

#add RNA mean /RNA cv
RNA_metainfo=pd.read_csv("../../../output/01-youth/02-correlation_calculation/RNA_metainfo_data.csv")
merged_df = pd.merge(df, RNA_metainfo, on='gene_id') 
#add gene length/cpg number
gene_metainfo=pd.read_csv("../../../output/01-youth/02-correlation_calculation/gene_metainfo.csv")
merged_df2 = pd.merge(merged_df, gene_metainfo, on=['gene_id','var_dim']) 
# Reorder
cor=[f'cor{j}' for j in range(1, 5001)]
new_columns_order = ['gene_id', 'gene_name','Gene_length', 'Cpg_number','datatype', 'mc_type', 'var_dim','RNA.NA.ratio','RNA.mean', 'RNA.cv','paired_cell_num', 'Correlation', 'P.value', 'Adjusted.P.value','new.P.value', 'new.P.adjust']+cor
merged_df2 = merged_df2[new_columns_order]

merged_df2.to_csv('../../../output/01-youth/02-correlation_calculation/shuffled_RNA_DNA_correlation_result_all.csv', index=False)               


