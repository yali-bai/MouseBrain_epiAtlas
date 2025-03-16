#########    All "our" in the following code refers to Joint Cabernet.
import sys
import pandas as pd
import numpy as np
import anndata as ad
from intervaltree import IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from matplotlib.font_manager import fontManager, FontProperties

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

modification, state, region, RNA_origin=sys.argv[1:5]
# modification: 5mCG, 5hmCG
# state: hyper, hypo
# region: promoter, genebody
# RNA_origin: our, zeng

path = f"{indir}/arial.ttf"
fontManager.addfont(path)

prop = FontProperties(fname=path)
sns.reset_defaults()
sns.set(context=None,
    style=None,
    palette=None,
    font=prop.get_name(),
    font_scale=None,
    color_codes=None,
    rc=None,)


with open("../../../../input/01-youth/subclass_order_for_integration_with_zeng.txt", "rt") as f:
    subclass_list=f.read().split("\n")[:-1]
subclass_list=[i for i in subclass_list if "NN" not in i]

df_RNA=pd.read_csv(f"{indir}/02-{RNA_origin}_subclass_mean_dat_final.csv",index_col=0).T

max_threshold_=0.2
cv_threshold_=0.5
df_RNA_selected=df_RNA.loc[subclass_list,:]
df_RNA_selected=df_RNA_selected.loc[subclass_list,
    (df_RNA_selected.max(axis=0)>max_threshold_) &
    (df_RNA_selected.apply(lambda col: col.std() / col.mean(), axis=0)>cv_threshold_)
            ]


DMR_states_df=ad.read_h5ad(f"{indir}/{modification}_DMR_states_neuron_1vsOthers.h5ad").to_df()


def get_gene_id(region_, trees, overlap_threshold=80):
    chrom, start_, end_ = region_.split('_')
    start_, end_ = int(start_), int(end_)
    intervals = trees[chrom][start_:end_]
    gene_id_list=[]
    for interval in intervals:
        a,b,c,d=start_, interval.begin, end_, interval.end
        if (min(c,d)-max(a,b))>=max(0.6*min((c-a),(d-b)),overlap_threshold):
            gene_id_list.append((min(c,d)-max(a,b),interval.data))
    if gene_id_list:
        return sorted(gene_id_list, key=lambda x : x[0], reverse=True)[0][1]
    else:
        return None


if region == "promoter":
    genes_df = pd.read_csv('../../../../input/reference_genome/PromoterUp2k.mm10.bed', sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id', 'gene_name', 'strand'])
else:
    genes_df = pd.read_csv('../../../../input/reference_genome/Genebody.mm10.bed', sep='\t', header=None, names=['chromosome', 'start', 'end', 'gene_id', 'gene_name', 'strand'])

df_frac_segment_subclass=dict()
df_frac_segment_subclass["5hmCG"]=ad.read_h5ad(f"{indir}/5hmC_merged_allc_with_{modification[:-1]}_segment.h5ad").to_df().iloc[:len(subclass_list),:]
df_frac_segment_subclass["total5mCG"]=ad.read_h5ad(f"{indir}/5mC_merged_allc_with_{modification[:-1]}_segment.h5ad").to_df().iloc[:len(subclass_list),:]
df_frac_segment_subclass["5hmCG"].index=df_RNA_selected.index
df_frac_segment_subclass["total5mCG"].index=df_RNA_selected.index
df_frac_segment_subclass["true5mCG"]=df_frac_segment_subclass["total5mCG"]-df_frac_segment_subclass["5hmCG"]
df_frac_segment_subclass["true5mCG"].clip(lower=0, inplace=True)

trees = {}
df_RNA_temp=df_RNA_selected.copy()
genes_df_temp=genes_df.loc[genes_df["gene_id"].isin(df_RNA_temp.columns)]
for chrom in ["chr"+str(i) for i in range(1,20,1)]:
    trees[chrom] = IntervalTree()
    sub_df = genes_df_temp[genes_df_temp['chromosome'] == chrom]
    for row in sub_df.iterrows():
        trees[chrom][row[1]["start"]:row[1]["end"]] = row[1]["gene_id"]

if state == "hyper":
    flag=1
elif state == "hypo":
    flag=-1


if region=="genebody":
    overlap_threshold_=500
elif region=="promoter":
    overlap_threshold_=100


temp_list=[]
for region_ in DMR_states_df.columns:
    diff_subclass=DMR_states_df.index[DMR_states_df[region_]==flag]
    gene_id=get_gene_id(region_,trees, overlap_threshold=overlap_threshold_)
    if diff_subclass.any() and gene_id:
        temp_list.append([region_,diff_subclass,gene_id])
df_DMR_geneid_diff=pd.DataFrame(
    sorted(
        temp_list,
        key=lambda x : sorted([subclass_list.index(i) for i in x[1]])
    ),
    columns=["segment", "diff_subclass", "gene_id"]
    )

df_DMR_geneid_diff=df_DMR_geneid_diff.groupby('gene_id').head(30)
# This leaves the first 30 DMRs for every gene.

df_DMR_geneid_diff["chromosome"]=df_DMR_geneid_diff["segment"].apply(lambda i : i.split("_")[0])
df_DMR_geneid_diff["start"]=df_DMR_geneid_diff["segment"].apply(lambda i : i.split("_")[1])
df_DMR_geneid_diff["end"]=df_DMR_geneid_diff["segment"].apply(lambda i : i.split("_")[2])

df_DMR_geneid_diff["diff_subclass_index"]=df_DMR_geneid_diff["diff_subclass"].apply(lambda x : sorted([subclass_list.index(i) for i in x]) )

df_DMR_geneid_diff["diff_subclass"]=df_DMR_geneid_diff["diff_subclass"].apply(lambda i : list(i))

df_DMR_geneid_diff.to_csv(f"{outdir}/{"_".join([modification, state, region, RNA_origin])}-DMR_GeneID_list.csv", index=False)

# with open(f"./DMR_GeneID_list/{"_".join([modification, state, region, RNA_origin])}-DMR_GeneID_list.txt", "wt") as f:
#     for i in df_DMR_geneid_diff['gene_id'].unique():
#         f.write(i+"\n")

df_frac_selected=dict()
df_frac_selected["5hmCG"]=df_frac_segment_subclass["5hmCG"].loc[:, df_DMR_geneid_diff['segment']]
df_frac_selected["total5mCG"]=df_frac_segment_subclass["total5mCG"].loc[:, df_DMR_geneid_diff['segment']]
df_frac_selected["true5mCG"]=df_frac_segment_subclass["true5mCG"].loc[:, df_DMR_geneid_diff['segment']]
temp1=df_RNA_temp.loc[:, df_DMR_geneid_diff['gene_id']]
temp2=temp1.div(temp1.mean(axis=0))
df_frac_selected["5hmCG"].index=temp2.index
df_frac_selected["total5mCG"].index=temp2.index
df_frac_selected["true5mCG"].index=temp2.index


font_size=5


for modification_show in ["5hmCG", "total5mCG", "true5mCG"]:


    plt.rcParams['figure.constrained_layout.use'] = True
    plt.figure(figsize=(2.5/2.54,6.0/2.54))
    ax=sns.heatmap(df_frac_selected[modification_show].T,vmin=(0 if modification_show != "total5mCG" else 0.5),vmax=(0.5 if modification_show == "5hmCG" else 1), cmap="rocket", yticklabels=False,xticklabels=False,
                  cbar_kws={"drawedges": False, "pad":0.01, 'shrink': 0.22, 'location': 'bottom',"orientation": "horizontal", "aspect": 6})
    ax.tick_params(labelsize=5,pad=1,length=0)
    cbar=ax.collections[0].colorbar
    cbar.ax.set_xticks([(0 if modification_show != "total5mCG" else 0.5),(0.5 if modification_show == "5hmCG" else 1)])
    cbar.ax.tick_params(labelsize=5,labelcolor="black",length=1,pad=1)
    
    plt.xlabel("neuronal subclasses",fontdict={"size":font_size}, labelpad=2)
    plt.ylabel(f"{state}-{"DHMR" if modification=="5hmCG" else "DMR"}s overlaping {"genebodies" if region =="genebody" else "promoters"}",fontdict={"size":font_size}, labelpad=2)
    plt.title(f"{"5hmCG" if modification_show=="5hmCG" else ("5mCG+5hmCG" if modification_show=="total5mCG" else "5mCG")} levels",fontdict={"size":font_size},pad=3)
    plt.savefig(f"{outdir}/DMR_plot-{"_".join([modification, state, region, RNA_origin])}-{modification_show}_levels.pdf")
    
    

    plt.rcParams['figure.constrained_layout.use'] = True
    plt.figure(figsize=(2.5/2.54,6.0/2.54))
    ax=sns.heatmap(df_frac_selected[modification_show].T,vmin=(0 if modification_show != "total5mCG" else 0.5),vmax=(0.5 if modification_show == "5hmCG" else 1), cmap="rocket", yticklabels=False,xticklabels=False,
                  cbar_kws={"drawedges": False, "pad":0.01, 'shrink': 0.22, 'location': 'bottom',"orientation": "horizontal", "aspect": 6})
    ax.tick_params(labelsize=5,pad=1,length=0)
    cbar=ax.collections[0].colorbar
    cbar.ax.set_xticks([(0 if modification_show != "total5mCG" else 0.5),(0.5 if modification_show == "5hmCG" else 1)])
    cbar.ax.tick_params(labelsize=5,labelcolor="black",length=1,pad=1)
    
    plt.xlabel("neuronal subclasses",fontdict={"size":font_size}, labelpad=2)
    plt.ylabel(f"{state}-{"DHMR" if modification=="5hmCG" else "DMR"}s overlaping {"genebodies" if region =="genebody" else "promoters"}",fontdict={"size":font_size}, labelpad=2)
    plt.title(f"{"5hmCG" if modification_show=="5hmCG" else ("5mCG+5hmCG" if modification_show=="total5mCG" else "5mCG")} levels",fontdict={"size":font_size},pad=3)
    plt.savefig(f"{outdir}/DMR_plot-{"_".join([modification, state, region, RNA_origin])}-{modification_show}_levels.png", dpi=600, bbox_inches="tight")





plt.rcParams['figure.constrained_layout.use'] = True
plt.figure(figsize=(2.5/2.54,6.0/2.54))
ax=sns.heatmap(temp2.T,vmin=0.2,vmax=1.8,cmap="cividis",yticklabels=False,xticklabels=False,
              cbar_kws={"drawedges": False, "pad":0.01, 'shrink': 0.22, 'location': 'bottom',"orientation": "horizontal", "aspect": 6})
ax.tick_params(labelsize=5,pad=1,length=0)
cbar=ax.collections[0].colorbar
cbar.ax.set_xticks([0.2,1.8])
cbar.ax.tick_params(labelsize=5,labelcolor="black",length=1,pad=1)
plt.xlabel("neuronal subclasses",fontdict={"size":font_size}, labelpad=2)
plt.ylabel(f"{state}-{"DHMR" if modification=="5hmCG" else "DMR"}s overlaping {"genebodies" if region =="genebody" else "promoters"}",fontdict={"size":font_size}, labelpad=2)
plt.title("normalized RNA expression",fontdict={"size":font_size},pad=3)
plt.savefig(f"{outdir}/RNA_plot_{"_".join([modification, state, region, RNA_origin])}.pdf")


plt.rcParams['figure.constrained_layout.use'] = True
plt.figure(figsize=(2.5/2.54,6.0/2.54))
ax=sns.heatmap(temp2.T,vmin=0.2,vmax=1.8,cmap="cividis",yticklabels=False,xticklabels=False,
              cbar_kws={"drawedges": False, "pad":0.01, 'shrink': 0.22, 'location': 'bottom',"orientation": "horizontal", "aspect": 6})
ax.tick_params(labelsize=5,pad=1,length=0)
cbar=ax.collections[0].colorbar
cbar.ax.set_xticks([0.2,1.8])
cbar.ax.tick_params(labelsize=5,labelcolor="black",length=1,pad=1)
plt.xlabel("neuronal subclasses",fontdict={"size":font_size}, labelpad=2)
plt.ylabel(f"{state}-{"DHMR" if modification=="5hmCG" else "DMR"}s overlaping {"genebodies" if region =="genebody" else "promoters"}",fontdict={"size":font_size}, labelpad=2)
plt.title("normalized RNA expression",fontdict={"size":font_size},pad=3)
plt.savefig(f"{outdir}/RNA_plot_{"_".join([modification, state, region, RNA_origin])}.png",dpi=600, bbox_inches="tight")


print(f"number of DMR regions:\t{df_frac_selected["5hmCG"].shape[1]}")
print(f"number of corresponding genes:\t{len(temp2.T.index.unique())}")






