#!/usr/bin/env python
# coding: utf-8

##### 01.import packages #####
import anndata
import pandas as pd
from scipy import stats
import math
from scipy.stats import mannwhitneyu
import sys
import os
import pathlib
import pandas as pd
import seaborn as sns
from ALLCools.mcds import MCDS
from ALLCools.dataset import ALLCoolsDataset
import os
import matplotlib.pyplot as plt
import sys

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
# os.chdir("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/05.significant_DMR_DHMR_mcds/20250217_update.filter_before_calculating_p_value_adjust")

##### 03.read cell metainfo in #####
celldata = pd.read_csv("../../output/03-aging/mouse_young_and_old.metainfo.csv",low_memory=False,index_col=0)

##### 04.extract methylation levels #####
## 5mCG_5hmCG ##
mcds_path_list = list(pathlib.Path(indir).glob('TSO-joint-RNA_Mouse_5mC_all_cells.*.mcds'))

## set cell_name ##
celldata['mC_new'] = celldata['mC']
celldata.loc[celldata['age']=='old','mC_new'] += '.mm10.dna.tsv.gz'

## QC ##
mcds = MCDS.open(mcds_path_list, var_dim='segment', use_obs='allc_'+celldata[(celldata['total_QC'] == 1)]['mC_new'])

## generate matrix ##
var_dim='segment'
mcds = mcds.remove_chromosome(exclude_chromosome=["chrX", "chrY", "chrM", "chrL"], var_dim=var_dim)
mcds[var_dim + '_da_frac'] = mcds[var_dim + '_da'].sel(count_type="mc")/mcds[var_dim + '_da'].sel(count_type="cov")
df = pd.DataFrame(mcds[f'{var_dim}_da_frac'][:,:].values.reshape((mcds.dims['cell'],mcds.dims[var_dim])), index=mcds[f'{var_dim}_da_frac'].cell,columns=mcds[f'{var_dim}_da_frac'].segment)

## extract subclass and age information ##
df['mC_new']=df.index.str.replace("allc_", "")
merged_df = df.merge(celldata[celldata['total_QC'] == 1][['mC_new', 'lt_twice_subclass','age']], on='mC_new', how='left')
merged_df 

## create unique id for calculating mean value ##
merged_df['unique'] = merged_df['lt_twice_subclass']+'_'+merged_df['age']

## calculate mean value per subclass per age ##
mC_subclass_age_rowmeans = merged_df.groupby(['unique']).mean()
mC_subclass_age_rowmeans
mC_subclass_age_rowmeans['unique'] = mC_subclass_age_rowmeans.index

## deformation data ##
df_melted = mC_subclass_age_rowmeans.melt(id_vars='unique', var_name="segment", value_name='value')
df_melted

## save result ##
df_melted.to_csv(f"{outdir}/aging_DMR_DHMR_5mCG_5hmCG_mean_methy_level_of_subclass_diff_age.csv",index=True)

## 5hmCG ##
mcds_path_list = list(pathlib.Path(indir).glob('TSO-joint-RNA_Mouse_5hmC_all_cells.*.mcds'))

## set cell_name ##
celldata['hmC_new'] = celldata['hmC']
celldata.loc[celldata['age']=='old','hmC_new'] += '.mm10.dna.tsv.gz'

## QC ##
hmC_mcds = MCDS.open(mcds_path_list, var_dim='segment', use_obs='allc_'+celldata[(celldata['total_QC'] == 1)]['hmC_new'])

## ## generate matrix ##
var_dim='segment'
hmC_mcds = hmC_mcds.remove_chromosome(exclude_chromosome=["chrX", "chrY", "chrM", "chrL"], var_dim=var_dim)
hmC_mcds[var_dim + '_da_frac'] = hmC_mcds[var_dim + '_da'].sel(count_type="mc")/hmC_mcds[var_dim + '_da'].sel(count_type="cov")
hmC_df = pd.DataFrame(hmC_mcds[f'{var_dim}_da_frac'][:,:].values.reshape((hmC_mcds.dims['cell'],hmC_mcds.dims[var_dim])), index=hmC_mcds[f'{var_dim}_da_frac'].cell,columns=hmC_mcds[f'{var_dim}_da_frac'].segment)

## extract subclass and age information ##
hmC_df['hmC_new']=hmC_df.index.str.replace("allc_", "")
hmC_merged_df = hmC_df.merge(celldata[celldata['total_QC'] == 1][['hmC_new', 'lt_twice_subclass','age']], on='hmC_new', how='left')
hmC_merged_df 

## create unique id for calculating mean value ##
hmC_merged_df['unique'] = hmC_merged_df['lt_twice_subclass']+'_'+hmC_merged_df['age']

## calculate mean value per subclass per age ##
hmC_subclass_age_rowmeans = hmC_merged_df.groupby(['unique']).mean()
hmC_subclass_age_rowmeans
hmC_subclass_age_rowmeans['unique'] = hmC_subclass_age_rowmeans.index

## deformation data ##
hmC_df_melted = hmC_subclass_age_rowmeans.melt(id_vars='unique', var_name="segment", value_name='value')

## save result ##
hmC_df_melted.to_csv(f"{outdir}/aging_DMR_DHMR_5hmCG_mean_methy_level_of_subclass_diff_age.csv",index=True)

## true mC ##
## sort dataframe for calculating true mCG ##
df_reindex = df.reindex('allc_'+celldata['mC_new'])
hmC_df_reindex = hmC_df.reindex('allc_'+celldata['hmC_new'])
columns = list(df_reindex.columns)
columns[-1] = "hmC_new"
df_reindex.columns = columns

## test if rownames and colnames are matched ##
index_match = df_reindex.index.equals(hmC_df_reindex.index)
index_match  ### True

column_match = df_reindex.columns.equals(hmC_df_reindex.columns)
column_match  ### True

## The rownames and colnames of two datafrmaes should be the same when two dataframes are subtracted ##
hmC_df_reindex.index = celldata['uniq_id']
df_reindex.index = celldata['uniq_id']

## remove non-numeric column ##
df_reindex_backup = df_reindex.drop(columns=["hmC_new"])
hmC_df_reindex_backup = hmC_df_reindex.drop(columns=["hmC_new"])
#del df_reindex
#del hmC_df_reindex
true_mC_df = df_reindex_backup - hmC_df_reindex_backup

## extract subclass and age information ##
celldata['unique'] = celldata['uniq_id']
true_mC_df['unique']=true_mC_df.index
true_mC_merged_df = true_mC_df.merge(celldata[['unique', 'lt_twice_subclass','age']], on='unique', how='left')
true_mC_merged_df 
true_mC_merged_df['Unique_ID'] =  true_mC_merged_df['unique']

## create unique id for calculating mean value ##
true_mC_merged_df['unique'] = true_mC_merged_df['lt_twice_subclass']+'_'+true_mC_merged_df['age']

## calculate mean value per subclass per age ##
true_mC_subclass_age_rowmeans = true_mC_merged_df.groupby(['unique']).mean()

## deformation data ##
true_mC_subclass_age_rowmeans['unique'] = true_mC_subclass_age_rowmeans.index
true_mC_df_melted = true_mC_subclass_age_rowmeans.melt(id_vars='unique', var_name="segment", value_name='value')

## save result ##
true_mC_df_melted.to_csv(f"{outdir}/aging_DMR_DHMR_true_5mCG_mean_methy_level_of_subclass_diff_age.csv",index=True)

