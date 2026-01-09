################02Feature Basic Filtering#########################
##################################################################
##################################################################
import pathlib
import pandas as pd
import seaborn as sns
from ALLCools.mcds import MCDS
from ALLCools.dataset import ALLCoolsDataset
import os
import matplotlib.pyplot as plt
import sys
import numpy as np
import anndata

dt = sys.argv[1]
os.chdir("./")

##### mc QC
celldata = pd.read_csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",low_memory=False)
if dt == "mC":
    colname_qc = 'mC_QC'
    colname_name = 'mC_SampleID'
else:
    colname_qc = 'hmC_QC'
    colname_name = 'hmC_SampleID'

QC_cells = ('allc_'+celldata[(celldata[colname_qc] == 1)][colname_name])
var_dim = 'segment'
mcds_path_list = list(pathlib.Path('./').glob('TSO-joint-RNA_Mouse_'+dt+'_all_cells.*.mcds'))
mcds = MCDS.open(mcds_path_list, var_dim=var_dim, use_obs=QC_cells)
mcds
mcds = mcds.remove_chromosome(exclude_chromosome=["chrX", "chrY", "chrM", "chrL"], var_dim=var_dim)
mcds[var_dim + '_da_frac'] = mcds[var_dim + '_da'].sel(count_type="mc")/mcds[var_dim + '_da'].sel(count_type="cov")

cg_adata = mcds.get_adata(mc_type="CGN",
                                    var_dim = var_dim,
                                    select_hvf=False)

df_5mCG=cg_adata.to_df().T
for chr_number in range(1, 20, 1): 
    temp_df=df_5mCG[df_5mCG.index.str.startswith(f"chr{chr_number}_")]
    temp_adata = anndata.AnnData(temp_df.T)
    temp_adata.write(f'{dt}G_frac_segment_by_cell_chr{chr_number}.h5ad')



