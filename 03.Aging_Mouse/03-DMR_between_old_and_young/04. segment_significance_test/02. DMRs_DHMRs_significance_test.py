##### 01.import packages #####
import anndata  
import pandas as pd
from scipy import stats 
import math  
from scipy.stats import mannwhitneyu  
import sys
import os
from scipy.stats import kstest
from scipy.stats import wilcoxon
import numpy as np

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
os.chdir("./")

##### 03.define function of permutation_test #####
def permutation_test(data1, data2, num_permutations=1000, statistic=np.mean):
    # Calculate the difference of the original statistics
    obs_diff = np.abs(statistic(data1) - statistic(data2))
    
    # define Statistical differences after storage replacement
    perm_diffs = []
    
    # Perform a replacement test
    for _ in range(num_permutations):
        # Randomly replace data labels
        combined = np.concatenate([data1, data2])
        np.random.shuffle(combined)
        
        # Redivide the dataset
        perm_data1 = combined[:len(data1)]
        perm_data2 = combined[len(data1):]
        
        # Calculate the difference of statistics after replacement
        perm_diff = np.abs(statistic(perm_data1) - statistic(perm_data2))
        perm_diffs.append(perm_diff)
    
    # calculate p-value 
    p_value = np.sum(np.array(perm_diffs) >= obs_diff) / num_permutations
    
    return obs_diff, p_value

##### 03.define function of significance test #####
def DMR_diff_compute(datatype,chr_number,n1,n2):
    info=pd.read_csv("../../../04.data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv",low_memory=False)
    info['subclass'] = info['subclass_label']
    mc=anndata.read_h5ad(f"{indir}/{datatype}_frac_segment_by_cell_chr{chr_number}.h5ad").to_df()
    info=info[info["total_QC"]==1]
    if datatype == "5mCG":
        info.index='allc_'+info['mC_ID']
    else:
        info.index='allc_'+info['hmC_ID']

    gene_results_all = pd.DataFrame()
    n1=int(n1)
    n2=int(n2)
    for nrow in range(n1-1,n2):   #0-base
        for cl in info['subclass'].unique():
            print(f"{nrow}-{cl}")
            y_cells=info.index[(info['old_young'] == 'young') & (info['subclass'] == cl)] # index of cells with candidate subclass and young age
            o_cells=info.index[(info['old_young'] == 'old') & (info['subclass'] == cl)] # index of cells with candidate subclass and old age
            gene=mc.columns[nrow] # candidate segment
            y_values = mc.loc[mc.index.isin(y_cells),gene] # young cellnames 
            o_values = mc.loc[mc.index.isin(o_cells +".mm10.dna.tsv.gz"),gene] # old cellnames 
            y_na_ratio=(y_values.isna().sum())/(y_values.shape[0]) # non na cell ratio in young samples
            o_na_ratio=(o_values.isna().sum())/(o_values.shape[0]) # non na cell ratio in old samples
            y_clean=y_values.dropna().values # methylation levels in young samples after removing na
            o_clean=o_values.dropna().values # methylation levels in old samples after removing na
            if ((y_clean!=0).sum()>2)&((o_clean!=0).sum()>2):
                _, ttest_p = stats.ttest_ind(y_clean,o_clean,equal_var=False,nan_policy="omit") # t test
                _, mannwhitneyu_p = mannwhitneyu(y_clean,o_clean, alternative='two-sided') # wilcox test
                _, ks_p =stats.ks_2samp(y_clean,o_clean) # ks test
                _, permutation_p = permutation_test(y_clean,o_clean) # permutation test
                y_mean=y_clean.mean() # mean methylation levels in young samples after removing na
                o_mean=o_clean.mean() # mean methylation levels in old samples after removing na
                y_num = len(y_clean) # non na cell number in young samples
                o_num = len(o_clean) # non na cell number in old samples
                logFC=math.log2(o_mean/y_mean) # logFC
                df=pd.DataFrame({
                'chrom':gene,
                'level':"subclass",
                'cluster':cl,
                'young_mean':y_mean,
                'old_mean':o_mean,
                'young_number':y_num,
                'old_number':o_num,
                'young_na_ratio':y_na_ratio,
                'old_na_ratio':o_na_ratio,
                'logFC':logFC,
                't.test_p':ttest_p,
                'mannwhitneyu_p':mannwhitneyu_p,
                'ks_p':ks_p,
                'permutation_p':permutation_p
                },
                index=[0])
                gene_results_all = pd.concat([gene_results_all, df], axis=0, ignore_index=True) 
    gene_results_all.to_csv(f"{outdir}/{datatype}_DMR_DHMR.chr{chr_number}_{n1}_{n2}.csv",index=False)


if __name__ == "__main__":
    args = sys.argv[1:]
    print(args)
    
    datatype = args[0]
    chr_number = args[1]
    n1 = args[2]
    n2 = args[3]
    ## run function ##
    DMR_diff_compute(datatype,chr_number,n1,n2)
