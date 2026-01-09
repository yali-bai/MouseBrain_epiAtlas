import anndata
import pandas as pd
import numpy as np
import pathlib
from ALLCools.mcds import RegionDS, MCDS
import xarray as xr
import yaml
import subprocess
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.font_manager import fontManager, FontProperties

# "indir" is a custom input path, and "outdir" is a custom output path.
indir = './RegionDS_data'
outdir = './plot'

path = "/arial/arial.ttf"
fontManager.addfont(path)
plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.size'] = 5

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


list_RegionDS=list(pathlib.Path(indir).glob("*DMR_states*"))

print(list_RegionDS)


for i in list_RegionDS:
    for state_ in ["hyper", "hypo", "all"]:
        print(state_, i)
        dmr_ds = RegionDS.open(i)
        array_dmr_state=dmr_ds["dmr_state"].values
        array_annotation=dmr_ds["dmr_genome-features_da"].values
        if state_ == "all":
            array_dmr_annotation=array_annotation[(np.abs(array_dmr_state)==1).sum(axis=1)>0]
        elif state_ == "hyper":
            array_dmr_annotation=array_annotation[(array_dmr_state==1).sum(axis=1)>0]
        elif state_ == "hypo":
            array_dmr_annotation=array_annotation[(array_dmr_state==-1).sum(axis=1)>0]
        else:
            raise Exception("! wrong state_ !")
        array_dmr_annotation=array_dmr_annotation[array_dmr_annotation.sum(axis=1)>0]
        array_dmr_annotation_excluded=array_dmr_annotation.copy()
        array_dmr_annotation_excluded[array_dmr_annotation[:,2:].sum(axis=1)>0,:2]=0
        array_dmr_annotation_excluded_normalized=array_dmr_annotation_excluded/array_dmr_annotation_excluded.sum(axis=1,keepdims=True)
        annotation_Series=pd.Series(array_dmr_annotation_excluded_normalized.sum(axis=0),
             index=dmr_ds["genome-features"].values)
        annotation_Series_percentage=annotation_Series*100/annotation_Series.sum()

        # Data
        labels=annotation_Series_percentage.index
        #colors = ['#4B0082', '#3F3B99', '#2EBAE4', '#8FE67A', '#51B96E', '#E5E059', '#FFDB3D']
        colors=["#3498DB", "#FFDC75", "#E74C3C", "#9B59B6", "#2ECC71", "#34495E"]

        # Create a pie chart
        plt.rcParams['figure.constrained_layout.use'] = False
        plt.figure(figsize=(53/25.4, 30/25.4))
        wedges, texts, autotexts = plt.pie(annotation_Series_percentage, autopct='',
                        startangle=90, colors=colors,radius=0.5,
                        wedgeprops={'linewidth': 0., 'edgecolor': 'white'})

        # Set the aspect ratio to be equal so the pie will be circular
        plt.axis('equal')

        # Add legend with labels and percentages
        plt.legend(wedges, [f'{label} ({size:.2f}%)' for label, size in zip(labels, annotation_Series_percentage)],
                   loc="center left",
                   bbox_to_anchor=(1, 0, 0.5, 1),
                  fontsize=5)

        plt.tight_layout(pad=0.1)
        FigName=state_+"-"+str(i).split("/")[-1].replace("_DMR_states_","-")+"-pie_chart.pdf"
        plt.savefig(f"{outdir}/{FigName}")
        plt.savefig(f"{outdir}/{FigName.replace('pdf', 'png')}",dpi=600)
        #plt.show()




















