# ----------------------------------------------------------------------------------------------------------------------------
# Script Name: 03-DNA_methylation_ratio_subclass_mean.py
# Author: yangf
# Created: 2025-10-23
# Description: Python script used to generate the DNA methylation ratio matrix of various combinations at the subclass level.
# ----------------------------------------------------------------------------------------------------------------------------

import sys
import pandas as pd
import numpy as np
import anndata as ad

def merge_mean(omics, element, group, subclass_labels):
    h5ad_file = f"5{omics}_{element}.{group}.pass_total_QC.h5ad"
    h5adObj = ad.read_h5ad(h5ad_file)
    counts = h5adObj.X
    counts = pd.DataFrame(counts, index=h5adObj.obs_names, columns=h5adObj.var_names)
    counts.index = counts.index.str.replace("allc_", "")
    
    csv_output_path = f"./01_3C/output/03-5{omics}_{element}.{group}.csv.gz"
    subclass_output_path = f"./01_3C/output/03-5{omics}_{group}_{element}_subclass_mean_dat_final.csv"
    
    counts.to_csv(csv_output_path, index=True, compression='gzip')
    
    metaInfo = pd.read_csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv")
    metaInfo = metaInfo[metaInfo['total_QC'] !=0]
    metaInfo_label = metaInfo[['subclass_label', f"{omics}_SampleID"]]
    
    counts_index = counts.index
    in_omics = counts_index.isin(metaInfo_label[f"{omics}_SampleID"])
    if all(in_omics):
        counts_filtered = counts
        print(f"All counts row names are in metaInfo_label['{in_omics}']. The original counts DataFrame is retained.")
    else:
        counts_filtered = counts[in_omics]
        print(f"Not all counts row names are in metaInfo_label['{in_omics}']. Only the rows that are in {in_omics} are retained.")
    
    counts_subclass = metaInfo_label.set_index(f"{omics}_SampleID").loc[counts_filtered.index, 'subclass_label']
    counts_filtered['subclass'] = counts_subclass
    counts_final = counts_filtered.groupby('subclass').agg(lambda x: np.nanmean(x, axis=0))
    counts_final = counts_final.T
    
    if all(col in subclass_labels for col in counts_final.columns):
        counts_final = counts_final[subclass_labels]
        print("All columns are in the subclass list. The columns have been sorted.")
    else:
        missing_columns = [col for col in counts_final.columns if col not in subclass_labels]
        print("Some columns are not in the subclass list:", missing_columns)
    
    counts_final.to_csv(subclass_output_path, index=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)

    omics = args[0]  # "hmC"
    element = args[1] # "chrom100k"
    group = args[2] # "CG"
    
    subclass = [
    'L2/3 IT CTX Glut','L4/5 IT CTX Glut','L5 IT CTX Glut','L6 IT CTX Glut','IT AON-TT-DP Glut','LA-BLA-BMA-PA Glut','L2/3 IT RSP Glut','L4 RSP-ACA Glut','L5 ET CTX Glut','SUB-ProS Glut','CA1-ProS Glut','CA3 Glut','CLA-EPd-CTX Car3 Glut','L5 NP CTX Glut','L6 CT CTX Glut','DG Glut','OB Eomes Ms4a15 Glut','OB-in Frmd7 Gaba','OB-out Frmd7 Gaba','Sncg Gaba','Lamp5 Gaba','Pvalb Gaba','Sst Gaba','STR D1 Gaba','STR D2 Gaba','ACB-BST-FS D1 Gaba','Astro-TE NN','Oligo NN','OPC NN','OEC NN','Microglia NN'
    ]

    merge_mean(omics, element, group, subclass_labels=subclass)
