# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

import sys
import pandas as pd
import numpy as np
import anndata as ad

def merge_mean(age,classs,omics, element, group):
    if age=="old":
        h5ad_file = f"{indir}/{omics}_{element}.{group}.pass_paired_QC.h5ad"
    else:
        h5ad_file = f"{indir}/{element}/{omics}_{element}.{group}.pass_paired_QC.h5ad"
    seuratObj = ad.read_h5ad(h5ad_file)
    counts = seuratObj.X
    counts = pd.DataFrame(counts, index=seuratObj.obs_names, columns=seuratObj.var_names)
    counts.index = counts.index.str.replace("allc_", "")
    metaInfo_file = f"../../04data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv"
    metaInfo = pd.read_csv(metaInfo_file,index_col=0)
    if omics=='5hmC':
        metaInfo_label = metaInfo[[f'{classs}_label', 'hmC_ID']][metaInfo['old_young']==age]
        counts_class = metaInfo_label.set_index('hmC_ID').loc[counts.index, f'{classs}_label']
    else:
        metaInfo_label = metaInfo[[f'{classs}_label', 'mC_ID']][metaInfo['old_young']==age]
        counts_class = metaInfo_label.set_index('mC_ID').loc[counts.index, f'{classs}_label']
    counts['class'] = counts_class
    counts_final = counts.groupby('class').agg(lambda x: np.nanmean(x, axis=0))
    counts_final = counts_final.T
    output_file = f"{outdir}/{classs}/{omics}_{group}_{element}_{classs}_mean_dat_final.csv"
    counts_final.to_csv(output_file)

if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)

    age = args[0]#"old"
    classs = args[1]#"subclass"
    omics = args[2]  # "5hmC"
    element = args[3] # "geneslop2k"
    group = args[4] # "CG"
    

    merge_mean(age,classs,omics, element, group)

# datatype=['5mC','5hmC','true_5mC']#
# elements=['genebody','promoter']
# groups=['CG','CH']
# classes=['class','subclass','three_class']

# for omics in datatype:
#     for element in elements:
#         for group in groups:
#             print(f"Running:{omics}_{element}_{group}")
#             merge_mean(omics, element, group)
