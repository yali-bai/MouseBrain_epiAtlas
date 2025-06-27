# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
import sys
import pandas as pd
import numpy as np
import anndata as ad
import rpy2.robjects as ro 
from rpy2.robjects import pandas2ri
pandas2ri.activate() 

def merge_mean(age,classs):
    rds_data = ro.r['readRDS'](f"{indir}/our_RNA_data_matrix.rds") 
    counts = pandas2ri.rpy2py(rds_data) 
    counts=counts.T
    counts.index = counts.index.str.replace(".*\\.\\._", "")
    metaInfo_file = f"../../04data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv"
    metaInfo = pd.read_csv(metaInfo_file,index_col=0)
    metaInfo_label = metaInfo[f'{classs}_label'][metaInfo['old_young']==age]
    age_counts=counts.loc[metaInfo_label.index,:]
    age_counts_class = metaInfo_label.loc[age_counts.index]
    age_counts['class'] = age_counts_class
    counts_final = age_counts.groupby('class').agg(lambda x: np.nanmean(x, axis=0))
    counts_final = counts_final.T
    output_file = f"{outdir}/{classs}/RNA_{classs}_mean_dat_final.csv"
    counts_final.to_csv(output_file)

if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)

    age = args[0]#"old"
    classs = args[1]#"subclass"


    merge_mean(age,classs)

# datatype=['5mC','5hmC','true_5mC']#
# elements=['genebody','promoter']
# groups=['CG','CH']
# classes=['class','subclass','three_class']

# for omics in datatype:
#     for element in elements:
#         for group in groups:
#             print(f"Running:{omics}_{element}_{group}")
#             merge_mean(omics, element, group)
