import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker

from matplotlib.font_manager import fontManager, FontProperties

path = "/arial/arial.ttf"
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

the_padding=0.6
red_alpha=0.7
hist_width=44
the_upper_ylim=1.075
bar_width=0.55
errorbar_setting={"ecolor":"black","elinewidth":0.6,"capsize":2.5,"capthick":0.6,"fmt":"none"}

data_5mC=pd.read_csv(r"mCH_data_for_5mC.txt",sep=" ",header=None,index_col=0)

data_5mC_hist=data_5mC.sum(axis=0)
data_5mC_hist[99]+=data_5mC_hist[100]
del data_5mC_hist[100]
data_5mC_hist_normalize2=data_5mC_hist/np.sum(data_5mC_hist)
data_5mC_hist_normalize2


plt.figure(figsize=(hist_width/25.4, 46/25.4))
ax=sns.barplot(x=[5*i for i in range(20)],y=np.sum(np.array(data_5mC_hist_normalize2).reshape((20,5)),axis=1),width=1,color = "#316e99")
ax.tick_params(labelsize=6)
plt.xticks([i-0.5 for i in range(0,21,5)],["0.0","0.25","0.5\nread mCH rate","0.75","1.0"],rotation='horizontal')
plt.yticks([5*i/100 for i in range(0,21,5)],rotation='horizontal')
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,19.5))
ax.set_title("Normalized read counts\nacross read mCH rates", fontdict={'size': 6})
ax.vlines([9.5], 0, the_upper_ylim, linestyles='dashed', colors='red',alpha=red_alpha)
plt.tight_layout(pad=the_padding)
plt.savefig(r"mCH_ratio_hist_5mC.pdf")
plt.show()


data_5hmC=pd.read_csv(r"mCH_data_for_5hmC.txt",sep=" ",header=None,index_col=0)

data_5hmC_hist=data_5hmC.sum(axis=0)
data_5hmC_hist[99]+=data_5hmC_hist[100]
del data_5hmC_hist[100]
data_5hmC_hist_normalize2=data_5hmC_hist/np.sum(data_5hmC_hist)
data_5hmC_hist_normalize2

plt.figure(figsize=(hist_width/25.4, 46/25.4))
ax=sns.barplot(x=[5*i for i in range(20)],y=np.sum(np.array(data_5hmC_hist_normalize2).reshape((20,5)),axis=1),width=1,color = "#316e99")
ax.tick_params(labelsize=6)
plt.xticks([i-0.5 for i in range(0,21,5)],["0.0","0.25","0.5\nread mCH rate","0.75","1.0"],rotation='horizontal')
plt.yticks([5*i/100 for i in range(0,21,5)],rotation='horizontal')
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,19.5))
ax.set_title("Normalized read counts\nacross read mCH rates", fontdict={'size': 6})
ax.vlines([9.5], 0, the_upper_ylim, linestyles='dashed', colors='red',alpha=red_alpha)
plt.tight_layout(pad=the_padding)
plt.savefig(r"mCH_ratio_hist_5hmC.pdf")
plt.show()

## 

batch_data=pd.read_csv(r"filter_QC.csv")
batch_data=batch_data[~batch_data["SampleID"].str.contains("384plate")]
batch_data_Ratios=batch_data.loc[:,["SampleID",'dna_reads_ratio','rna_reads_ratio']]
batch_data_Ratios.index=batch_data_Ratios["SampleID"]
batch_data_Ratios=batch_data_Ratios.loc[:,['dna_reads_ratio','rna_reads_ratio']]
def read_percentage(percentage_str):
    percentage_float = float(percentage_str.replace("%", "")) / 100.0
    return percentage_float
 
batch_data_Ratios_num=batch_data_Ratios.apply(lambda x : [read_percentage(i) for i in x])
batch_data_Ratios_num["Type"]=[ x.split("joint")[1].split("_")[0] for x in batch_data_Ratios_num.index]
batch_data_Ratios_num["dna_passing_filter_ratio"]=[dna/(dna+rna) for (dna,rna) in zip(batch_data_Ratios_num["dna_reads_ratio"],batch_data_Ratios_num["rna_reads_ratio"])]

plt.figure(figsize=((60-hist_width)/25.4, 46/25.4))
the_data=batch_data_Ratios_num[batch_data_Ratios_num["Type"].str.contains("5hmC")]
ax=sns.barplot(y="dna_passing_filter_ratio", width=bar_width, data=the_data, errorbar=None)
the_mean=np.mean(the_data["dna_passing_filter_ratio"])
the_std=np.std(the_data["dna_passing_filter_ratio"])
plt.errorbar(x=[0],y=the_mean,yerr=the_std,**errorbar_setting)
plt.xticks([0],["Reads ratio    \nafter CH filter    "],rotation='horizontal')
ax.tick_params(labelsize=6)
ax.set_ylabel("",fontdict={'size': '6'})
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,0.5))
ax.set_title(f"{the_mean:.2f} ± {the_std:.2f}    \n(mean ± SD)    ", fontdict={'size': 6})
plt.tight_layout(pad=the_padding)
plt.savefig(r"cell_5hmC_viable_ratio.pdf")
plt.show()


plt.figure(figsize=((60-hist_width)/25.4, 46/25.4))
the_data=batch_data_Ratios_num[batch_data_Ratios_num["Type"].str.contains("5mC")]
ax=sns.barplot(y="dna_passing_filter_ratio", width=bar_width, data=the_data, errorbar=None)
the_mean=np.mean(the_data["dna_passing_filter_ratio"])
the_std=np.std(the_data["dna_passing_filter_ratio"])
plt.errorbar(x=[0],y=the_mean,yerr=the_std,**errorbar_setting)
plt.xticks([0],["Reads ratio    \nafter CH filter    "],rotation='horizontal')
ax.tick_params(labelsize=6)
ax.set_ylabel("",fontdict={'size': '6'})
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,0.5))
ax.set_title(f"{the_mean:.2f} ± {the_std:.2f}    \n(mean ± SD)    ", fontdict={'size': 6})
plt.tight_layout(pad=the_padding)
plt.savefig(r"cell_5mC_viable_ratio.pdf")
plt.show()

##### RNA #####
data_RNA=pd.read_csv(r"mCH_data_for_RNA.txt",sep=" ",header=None,index_col=0)
data_RNA_hist=data_RNA.sum(axis=0)
data_RNA_hist[99]+=data_RNA_hist[100]
del data_RNA_hist[100]
data_RNA_hist_normalize2=data_RNA_hist/np.sum(data_RNA_hist)
data_RNA_hist_normalize2

plt.figure(figsize=(hist_width/25.4, 46/25.4))
ax=sns.barplot(x=[5*i for i in range(20)],y=np.sum(np.array(data_RNA_hist_normalize2).reshape((20,5)),axis=1),width=1,color = "#316e99")
ax.tick_params(labelsize=6)
plt.xticks([i-0.5 for i in range(0,21,5)],["0.0","0.25","0.5\nread mCH rate","0.75","1.0"],rotation='horizontal')
plt.yticks([5*i/100 for i in range(0,21,5)],rotation='horizontal')
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,19.5))
ax.set_title("Normalized read counts\nacross read mCH rates", fontdict={'size': 6})
ax.vlines([17.5], 0, the_upper_ylim, linestyles='dashed', colors='red',alpha=red_alpha)
plt.tight_layout(pad=the_padding)
plt.savefig(r"mCH_ratio_hist_RNA.pdf")
plt.show()

with open(r"cell_RNA_viable_ratio.txt", "rt") as f:
    batch_data_RNA_1=f.read()

batch_data_RNA_2=batch_data_RNA_1.split("\n")[:-1]
batch_data_RNA_3=[float(i) for i in batch_data_RNA_2]
len(batch_data_RNA_3)

#
plt.figure(figsize=((60-hist_width)/25.4, 46/25.4))

ax=sns.barplot(y=batch_data_RNA_3, width=bar_width, errorbar=None,color = "#316e99")
the_mean=np.mean(batch_data_RNA_3)
the_std=np.std(batch_data_RNA_3)
plt.errorbar(x=[0],y=the_mean,yerr=the_std,**errorbar_setting)
plt.xticks([0],["Reads ratio    \nafter CH filter    "],rotation='horizontal')

ax.tick_params(labelsize=6)
ax.set_ylabel("",fontdict={'size': '6'})
ax.set_ylim((0,the_upper_ylim))
ax.set_xlim((-0.5,0.5))
ax.set_title(f"{the_mean:.2f} ± {the_std:.2f}    \n(mean ± SD)    ", fontdict={'size': 6})
plt.tight_layout(pad=the_padding)
plt.savefig(r"cell_RNA_viable_ratio.pdf")
plt.show()