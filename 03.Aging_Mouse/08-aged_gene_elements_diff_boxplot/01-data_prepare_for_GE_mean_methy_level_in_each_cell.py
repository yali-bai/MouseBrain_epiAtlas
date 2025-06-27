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
import anndata as ad
import numpy as np

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 03.01. hmC #####
## young ##
young_promoter = ad.read_h5ad(f"{indir}/5hmC_promoter.CG.pass_paired_QC.h5ad")

matrix = young_promoter.X
cell_sum = np.nanmean(matrix,axis=1)
cell_sum 

df=pd.DataFrame(
    {'sampleid':young_promoter.obs_names,
    'gene_mean_level':cell_sum
    })
df

df.to_csv(f"{outdir}/young_promoter_hmCG.mean_methyl.csv")

## old ##
old_promoter = ad.read_h5ad(f"{indir}/5hmC_promoter.CG.pass_paired_QC.h5ad")
matrix = old_promoter.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':old_promoter.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/old_promoter_hmCG.mean_methyl.csv")

##### 03.02. mC #####
## young ##
young_promoter = ad.read_h5ad(f"{indir}/5mC_promoter.CG.pass_paired_QC.h5ad")
matrix = young_promoter.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':young_promoter.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/young_promoter_mCG.mean_methyl.csv")

## old ##
old_promoter = ad.read_h5ad(f"{indir}/5mC_promoter.CG.pass_paired_QC.h5ad")
matrix = old_promoter.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':old_promoter.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/old_promoter_mCG.mean_methyl.csv")

##### 04. genebody #####
##### 04.01. hmC #####
## young ##
young_genebody = ad.read_h5ad(f"{indir}/5hmC_genebody.CG.pass_paired_QC.h5ad")
matrix = young_genebody.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':young_genebody.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/young_genebody_hmCG.mean_methyl.csv")

## old ##
old_genebody = ad.read_h5ad(f"{indir}/5hmC_genebody.CG.pass_paired_QC.h5ad")
matrix = old_genebody.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':old_genebody.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/old_genebody_hmCG.mean_methyl.csv")

##### 04.02. mC #####
## young ##
young_genebody = ad.read_h5ad(f"{indir}/5mC_genebody.CG.pass_paired_QC.h5ad")
matrix = young_genebody.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':young_genebody.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv(f"{outdir}/young_genebody_mCG.mean_methyl.csv")

## old ##
old_genebody = ad.read_h5ad(f"{indir}/5mC_genebody.CG.pass_paired_QC.h5ad")
matrix = old_genebody.X
cell_sum = np.nanmean(matrix,axis=1)
df=pd.DataFrame(
    {'sampleid':old_genebody.obs_names,
    'gene_mean_level':cell_sum
    })
df.to_csv("old_genebody_mCG.mean_methyl.csv")

for GE in ['Enhancer','Intergenetic','UTR5','UTR3','Exon','Intro']:
    ##### 01. hmC #####
    ## young ##
    young_genebody = ad.read_h5ad(f"{indir}/5hmC_{GE}.CG.paired_QC.h5ad")
    matrix = young_genebody.X
    cell_sum = np.nanmean(matrix,axis=1)
    df=pd.DataFrame(
        {'sampleid':young_genebody.obs_names,
        'gene_mean_level':cell_sum
        })
    df.to_csv(GE+"_hmCG.mean_methyl.csv")
    
    
    ##### 02. mC #####
    ## young ##
    young_genebody = ad.read_h5ad(f"{indir}/5mC_{GE}.CG.paired_QC.h5ad")
    matrix = young_genebody.X
    cell_sum = np.nanmean(matrix,axis=1)
    df=pd.DataFrame(
        {'sampleid':young_genebody.obs_names,
        'gene_mean_level':cell_sum
        })
    df.to_csv(f"{outdir}/{GE}_mCG.mean_methyl.csv")
    

##### true mC #####
for GE in ['Enhancer','Intergenetic','UTR5','UTR3','Exon','Intro']:
    young_genebody_hmC = ad.read_h5ad(f"{indir}/5hmC_{GE}.CG.paired_QC.h5ad")
    matrix_hmC = young_genebody_hmC.X

    young_genebody_mC = ad.read_h5ad(f"{indir}/5mC_{GE}.CG.paired_QC.h5ad")
    matrix_mC = young_genebody_mC.X
    
    matrix_hmC_df = pd.DataFrame(matrix_hmC)
    matrix_mC_df = pd.DataFrame(matrix_mC)
    
    young_genebody_hmC.obs['cell'] = young_genebody_hmC.obs.index
    young_genebody_hmC.obs['sampleid'] = young_genebody_hmC.obs['cell'].apply(lambda x : x.replace('allc_', ''))
    young_genebody_hmC.obs['sampleid'] = young_genebody_hmC.obs['sampleid'].apply(lambda x : x.replace('.mm10.dna.tsv.gz', ''))
    young_genebody_hmC.obs
    
    metainfo = pd.read_csv("../../output/03.Aging_Mouse/mouse_young_and_old.metainfo.csv", index_col=0)
    metainfo.index = metainfo['hmC']
    matrix_hmC_df.index = metainfo.loc[young_genebody_hmC.obs['sampleid'],'uniq_id']
    matrix_hmC_df.index
    matrix_hmC_df.columns = young_genebody_hmC.var_names
    
    young_genebody_mC.obs['cell'] = young_genebody_mC.obs.index
    young_genebody_mC.obs['sampleid'] = young_genebody_mC.obs['cell'].apply(lambda x : x.replace('allc_', ''))
    young_genebody_mC.obs['sampleid'] = young_genebody_mC.obs['sampleid'].apply(lambda x : x.replace('.mm10.dna.tsv.gz', ''))
    young_genebody_mC.obs
    metainfo.index = metainfo['mC']
    matrix_mC_df.index = metainfo.loc[young_genebody_mC.obs['sampleid'],'uniq_id']
    matrix_mC_df = matrix_mC_df.reindex(matrix_hmC_df.index)
    matrix_mC_df.columns = young_genebody_mC.var_names
    
    if matrix_hmC_df.index.equals(matrix_mC_df.index) & matrix_hmC_df.columns.equals(matrix_mC_df.columns):
        matrix_true_mC_df = matrix_mC_df - matrix_hmC_df
        anndata_true_mC=ad.AnnData(matrix_true_mC_df.to_numpy())
        anndata_true_mC.var_names = matrix_mC_df.index
        anndata_true_mC.obs_names = matrix_mC_df.columns
        anndata_true_mC.write_h5ad(f"{outdir}/true_5mC_{GE}.CG.paired_QC.h5ad")

