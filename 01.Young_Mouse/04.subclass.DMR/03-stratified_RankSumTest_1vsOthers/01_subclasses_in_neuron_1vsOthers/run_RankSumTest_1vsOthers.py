import subprocess

for modification in ["5mC", "5hmC"]:
    for chr_number in range(1, 20, 1):
        print(modification,"chr",chr_number)
        if chr_number > 17:
            cpu=12
        elif chr_number > 7:
            cpu=8  
        else:
            cpu=6
            # Only 30 subclasses (vs Others), so only 30 processes needed.
            # However, the "others" cell-by-segment matrix is too big, probably about the size of  the whole matrix.
            # So we need to restrict the max_worker "cpu" parameter
        subprocess.run(f"sbatch -J neuron_RankSumTest_{modification}_chr{chr_number}_1vsOthers \
        -o ./logs/neuron_RankSumTest_{modification}_chr{chr_number}_1vsOthers.log \
        -e ./logs/neuron_RankSumTest_{modification}_chr{chr_number}_1vsOthers.log \
        --mem={"260G" if chr_number<10 else "240G"} \
        --partition=compute_fat \
        --cpus-per-task={cpu} \
        --time=150:00:00 \
           ./to_RankSumTest_1vsOthers.sh {modification} {chr_number} {cpu}",shell=True)
