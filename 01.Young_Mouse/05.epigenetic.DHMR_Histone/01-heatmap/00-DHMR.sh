import sys
import pandas as pd
import numpy as np
import anndata as ad

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

## 01-DHMR states
DMR_state = ad.read_h5ad(f"{indir}/5hmCG_DMR_states_neuron_1vsOthers.h5ad")
state_counts = DMR_state.X
state_counts = pd.DataFrame(state_counts, index=DMR_state.obs_names, columns=DMR_state.var_names)
na_values = state_counts.isna().sum().sum()
na_columns = state_counts.columns[state_counts.isna().all()]
state_counts.fillna(0, inplace=True)
#state_counts = state_counts.replace([0.0, -0.0, 1.0, -1.0], [0, 0, 1, -1])
state_counts = state_counts.applymap(lambda x: 0 if x == 0.0 or x == -0.0 else (1 if x == 1.0 else (-1 if x == -1.0 else x)))
state_counts_transposed = state_counts.T
#state_counts_transposed = state_counts_transposed.astype(int)
state_counts_transposed.to_csv(f'{outdir}/00-dhmr_state_collapsed_matrix_neuron.csv')

count_1_a = (state_counts_transposed == 1).sum()
count_neg1_a = (state_counts_transposed == -1).sum()
result_df = pd.DataFrame({
    'Subclass': state_counts_transposed.columns,
    'Hyper dhmr': count_1_a.values,
    'Hypo dhmr': count_neg1_a.values
})
result_df.to_csv(f'{outdir}/00-dhmr_hyper_hypo_subclass_number_neuron.csv')


## 02-specific_hyper_hypo_DHMR_statistics
import pandas as pd

df = pd.read_csv(f'{outdir}/00-dhmr_state_collapsed_matrix_neuron.csv', index_col=0)

# Create an empty dictionary to store the statistics.
specific_hyper_counts = {}
specific_hypo_counts = {}

for col in df.columns:
    # The current column is 1, and no other column is 1.
    specific_hyper = (df[col] == 1) & (df.drop(columns=col) != 1).all(axis=1)
    specific_hyper_counts[col] = specific_hyper.sum()
    # The current column is -1, and none of the other columns are -1.
    specific_hypo = (df[col] == -1) & (df.drop(columns=col) != -1).all(axis=1)
    specific_hypo_counts[col] = specific_hypo.sum()

specific_stats_df = pd.DataFrame({
    'Subclass': df.columns,
    'Specific_Hyper': list(specific_hyper_counts.values()),
    'Specific_Hypo': list(specific_hypo_counts.values())
})

specific_stats_df.to_csv(f'{outdir}/00-dhmr_specific_stats_df.csv')


## 03-DHMR-Genebody
import pandas as pd
from intervaltree import IntervalTree

df = pd.read_csv(f'{outdir}/00-dhmr_state_collapsed_matrix_neuron.csv', index_col=0)

genes_df = pd.read_csv('../../../input/reference_genome/Genebody.mm10.bed', sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id', 'gene_name', 'strand'])

# Create a dictionary to store the intervals of each chromosome.
trees = {}

for chrom in genes_df['chromosome'].unique():
    trees[chrom] = IntervalTree()
    sub_df = genes_df[genes_df['chromosome'] == chrom]
    for row in sub_df.itertuples():
        trees[chrom][row.start:row.end] = row.gene_name

# Define a function to find the interval and return the gene name.
def find_gene_name(region, trees):
    chrom, start, end = region.split('_')
    start, end = int(start), int(end)
    intervals = trees[chrom][start:end]
    for interval in intervals:
        if interval.begin <= start and interval.end >= end:
            return interval.data
    return None

# Add a new column to df.
df['gene_name'] = df.index.to_series().apply(lambda x: find_gene_name(x, trees))

# View the value count for the gene name column
gene_name_counts = df['gene_name'].value_counts(dropna=False)

# Extract the columns gene_name and region and remove the row where gene_name is None
df['region'] = df.index
gene_region_df = df[['gene_name', 'region']].copy()
gene_region_df = gene_region_df.dropna(subset=['gene_name'])

# Calculate the length of each region
gene_region_df['length'] = gene_region_df['region'].apply(lambda x: int(x.split('_')[2]) - int(x.split('_')[1]))

# Only the longest region is reserved for each gene name
longest_region_df = gene_region_df.loc[gene_region_df.groupby('gene_name')['length'].idxmax()]

longest_region_df = longest_region_df[['gene_name', 'region']]

# Reset index
longest_region_df = longest_region_df.reset_index(drop=True)

longest_region_df.to_csv(f'{outdir}/00-DHMR_in_genebody_paired.csv')

# Initializes the result dictionary
results = {
    'Subclass': [],
    'Specific_Hyper_Region': [],
    'Specific_Hyper_Count': [],
    'Specific_Hypo_Region': [],
    'Specific_Hypo_Count': []
}

for cell_type in df.columns[:-2]:  # Omit the gene name and region columns
    # Extract the Specific Hyper and Specific Hypo regions
    specific_hyper = df[(df[cell_type] == 1) & (~df.drop(columns=[cell_type, 'gene_name', 'region']).isin([1]).any(axis=1))][['region', 'gene_name']]
    specific_hypo = df[(df[cell_type] == -1) & (~df.drop(columns=[cell_type, 'gene_name', 'region']).isin([-1]).any(axis=1))][['region', 'gene_name']]
    results['Subclass'].append(cell_type)
    results['Specific_Hyper_Region'].append(specific_hyper['region'].tolist())
    results['Specific_Hyper_Count'].append(len(specific_hyper))
    results['Specific_Hypo_Region'].append(specific_hypo['region'].tolist())
    results['Specific_Hypo_Count'].append(len(specific_hypo))

specific_stats_df = pd.DataFrame(results)

specific_stats_df.to_csv(f'{outdir}/00-dhmr_specific_hyper_hypo_DHMR_include_region.csv')


##
results_2 = {
    'Subclass': [],
    'Hyper_Region': [],
    'Hyper_Count': [],
    'Hypo_Region': [],
    'Hypo_Count': []
}

for cell_type in df.columns[:-2]:  # Omit the gene name and region columns
    # Extract the Hyper and Hypo regions, other cell types are not considered
    hyper = df[df[cell_type] == 1][['region', 'gene_name']]
    hypo = df[df[cell_type] == -1][['region', 'gene_name']]
    
    # Append the results_2
    results_2['Subclass'].append(cell_type)
    results_2['Hyper_Region'].append(hyper['region'].tolist())
    results_2['Hyper_Count'].append(len(hyper))
    results_2['Hypo_Region'].append(hypo['region'].tolist())
    results_2['Hypo_Count'].append(len(hypo))

# Convert the results_2 into a DataFrame
Hyper_hypo_stats_df = pd.DataFrame(results_2)
Hyper_hypo_stats_df.to_csv(f'{outdir}/00-dhmr_hyper_hypo_DHMR_include_region.csv')
