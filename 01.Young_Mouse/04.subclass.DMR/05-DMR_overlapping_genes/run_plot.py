import subprocess

for modification in ["5mCG", "5hmCG"]:
    for state in ["hyper", "hypo"]:
        for region in ["promoter", "genebody"]:
            print(modification, state, region)
            subprocess.run(f"sbatch -J DMR_plot_{modification}_{state}_{region} \
        -o ./logs/DMR_plot_{modification}_{state}_{region}.log \
        -e ./logs/DMR_plot_{modification}_{state}_{region}.log \
        --mem=20G \
        --partition=compute_pro \
        --cpus-per-task=2 \
        --time=150:00:00 \
           ./do_plot_DMR_RNA.sh {modification} {state} {region}",shell=True)
