import subprocess

for modification in ["5mCG", "5hmCG"]:
    for state in ["hyper", "hypo"]:
        for region in ["promoter", "genebody"]:
            for RNA_origin in ["Joint_Cabernet", "zeng"]:
                print(modification, state, region, RNA_origin)
                subprocess.run(f"sbatch -J DMR_plot_{modification}_{state}_{region}_{RNA_origin} \
        -o ./logs/DMR_plot_{modification}_{state}_{region}_{RNA_origin}.log \
        -e ./logs/DMR_plot_{modification}_{state}_{region}_{RNA_origin}.log \
        --mem=20G \
        --partition=compute_pro \
        --cpus-per-task=2 \
        --time=150:00:00 \
           ./do_plot_DMR_RNA.sh {modification} {state} {region} {RNA_origin}",shell=True)
