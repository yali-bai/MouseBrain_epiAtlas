##### 01.import packages #####
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

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
os.chdir("./")

##### 03.5mCG_5hmCG #####
## read QC files in ##
old_celldata = pd.read_csv("../../../04.data/02.metainfo/03.Aging_Mouse/TSO-joint.mC_QC_stat.aged.csv",low_memory=False)

## QC and change sampleid ##
old_array = ('allc_'+old_celldata[(old_celldata['Library']=='mC') & (old_celldata['QC'] == 1) & (old_celldata['age'] == "old")]['SampleID']+'.mm10.dna.tsv.gz')
young_array = np.array('allc_'+old_celldata[(old_celldata['Library']=='mC') & (old_celldata['QC'] == 1) & (old_celldata['age'] == "young")]['SampleID'])
QC_cells = np.append(old_array, young_array)

var_dim = 'segment'
mcds_path_list = list(pathlib.Path(f'{indir}/03.segment_mcds/').glob('TSO-joint-RNA_Mouse_5mC_all_cells.*.mcds'))
mcds = MCDS.open(mcds_path_list, var_dim=var_dim, use_obs=QC_cells)
mcds

mcds = mcds.remove_chromosome(exclude_chromosome=["chrX", "chrY", "chrM", "chrL"], var_dim=var_dim)

## raw fraction is mc/cov ##
mcds[var_dim + '_da_frac'] = mcds[var_dim + '_da'].sel(count_type="mc")/mcds[var_dim + '_da'].sel(count_type="cov")

## convert to adata ##
cg_adata = mcds.get_adata(mc_type="CGN",
                                    var_dim = var_dim,
                                    select_hvf=False)

## output ##
#cg_adata.write_h5ad('subclass.5mC_' + var_dim + '.CG.pass_paired_QC.h5ad')

## output by chr ##
df_5mCG=cg_adata.to_df().T
for chr_number in range(1,20,1):
    temp_df=df_5mCG[df_5mCG.index.str.startswith(f"chr{chr_number}_")]
    temp_adata = anndata.AnnData(temp_df.T)
    temp_adata.write(f'{outdir}/5mCG_frac_segment_by_cell_chr{chr_number}.h5ad')


