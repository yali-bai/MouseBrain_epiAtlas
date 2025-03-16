import anndata
import pandas as pd
import numpy as np
import pathlib
from ALLCools.mcds import RegionDS, MCDS
import xarray as xr
import yaml
import subprocess
from tqdm import tqdm

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

def DMR_states_h5ad_to_region_ds(h5ad_path, output_dir):
    pathlib.Path(output_dir).mkdir(parents=True)
    with open(f"{output_dir}/.ALLCools", "w") as f:
        config = {"region_dim": "dmr", "ds_region_dim": {"dmr": "dmr"}}
        yaml.dump(config, f)
        
    df_DMR_states=anndata.read_h5ad(h5ad_path).to_df()
    df_DMR_states.index.name="sample"
    df_DMR_states.columns.name="dmr"
    
    dmr_ds=xr.Dataset({"dmr_state":df_DMR_states.T}, 
                coords={
                 "dmr_chrom": ("dmr",[i.split("_")[0] for i in df_DMR_states.columns]),
                "dmr_start": ("dmr",[i.split("_")[1] for i in df_DMR_states.columns]),
                "dmr_end": ("dmr",[i.split("_")[2] for i in df_DMR_states.columns]),
                 })
    
    dmr_ds.to_zarr(f"{output_dir}/dmr")
    subprocess.run(f"cp {indir}/chrom_sizes.txt {output_dir}/chrom_sizes.txt", shell=True)
    return


list_h5ad_path=list(pathlib.Path(f"{indir}/").glob(pattern="*/*DMR_states*.h5ad"))

list_output_path=[f"{outdir}/"+str(i).split("/")[-1].split(".")[0] for i in list_h5ad_path]

for q, w in tqdm(zip(list_h5ad_path, list_output_path)):
    print(q)
    print(w)
    print()

for q, w in tqdm(zip(list_h5ad_path, list_output_path), total=len(list_h5ad_path)):
    DMR_states_h5ad_to_region_ds(q, w)

list_RegionDS=list(pathlib.Path(outdir).glob("*"))

for i in tqdm(list_RegionDS):
    dmr_ds = RegionDS.open(i)
    dmr_ds.annotate_by_beds(slop=250,
                        bed_table=f'{indir}/new_genome_featue_bed.csv',
                        dim='genome-features',
                        bed_sorted=False,
                        cpu=38)



