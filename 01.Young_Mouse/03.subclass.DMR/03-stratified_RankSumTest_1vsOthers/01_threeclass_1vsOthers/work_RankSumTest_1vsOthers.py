import sys
import pandas as pd
import numpy as np
import anndata
from scipy import stats
from concurrent.futures import ProcessPoolExecutor, as_completed

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

modification, chr_number, cpu = sys.argv[1:4]
cpu=int(cpu)
modification2=modification[1:] # remove '5' from '5hmC' and '5mC'


def get_segment_length(segment_):
    chrom_, start_, end_ = segment_.split("_")
    return int(end_) - int(start_)

data=anndata.read_h5ad(f"{indir}/{modification}G_frac_segment_by_cell_chr{chr_number}.h5ad").to_df()
# now remove the segments whose length is shorter than 200bp
data_columns_1=pd.Series(data.columns)
data_columns_2=data_columns_1[data_columns_1.apply(get_segment_length)>=200]
data=data[data_columns_2]


meta_data=pd.read_csv("../../../../input/01-youth/RNA_DNA_match_name_QC_class_label.csv",header=0)
meta_data2=meta_data[meta_data["total_QC"]==1]
# meta_data3=meta_data2[["subclass_label", modification2]]
# meta_data_group=meta_data3.groupby("subclass_label")
meta_data3=meta_data2[["three_class_label", modification2]]
meta_data_group=meta_data3.groupby("three_class_label")


# with open("/share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/sample_info/04_order_for_class_subclass/subclass_order_for_integration_with_zeng.txt", "rt") as f:
#     subclass_list=f.read().split("\n")[:-1]


# #only keep NN
# subclass_list=[i for i in subclass_list if "NN" in i]

three_class_label_list=["Exc", "Inh", "Non"]

allc_name=lambda cellID : "allc_" + cellID

data_dict={three_class_ : data.loc[allc_name(meta_data_group.get_group(three_class_)[modification2])] for three_class_ in three_class_label_list}
data_columns=data.columns
del data # to save memory


# def get_other_subclass(subclass_):
#     other_subclass=[i for i in subclass_list if i != subclass_]
#     assert len(other_subclass)==len(subclass_list)-1
#     return pd.concat([data_dict[i] for i in other_subclass],axis=0)

def get_other_three_class(three_class_):
    other_three_class=[i for i in three_class_label_list if i != three_class_]
    assert len(other_three_class)==len(three_class_label_list)-1
    return pd.concat([data_dict[i] for i in other_three_class],axis=0)


def RankSumTest_1vsOthers(three_class_):
    result_df=pd.DataFrame(stats.mannwhitneyu(
                            data_dict[three_class_],
                            get_other_three_class(three_class_),
                            nan_policy="omit"
                            ))
    result_df.index=[f"{three_class_}_Others_statistic", f"{three_class_}_Others_pvalue"]
    return result_df

def main():
    total_result_list=[]
    with ProcessPoolExecutor(max_workers=cpu) as executor:
        for RankSumTest_result in executor.map(RankSumTest_1vsOthers, three_class_label_list):
            total_result_list.append(RankSumTest_result)
        total_result_df=pd.concat(total_result_list, axis=0)
        total_result_df.columns=data_columns
        temp_adata = anndata.AnnData(total_result_df)
        temp_adata.write(f'{outdir}/{modification}G_RankSumTest_result_chr{chr_number}.h5ad')


if __name__ == "__main__":
    main()
