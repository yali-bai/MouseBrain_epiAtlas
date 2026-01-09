import subprocess
import pandas as pd
import time

csv_path = "../../../03.data/02.metainfo/01.Young_Mouse/subclass_name.Joint_Cabernet_corresponding_to_3C.csv"
df = pd.read_csv(csv_path, sep=',')

for subclass_ in df["subclass_liu"]:
    print(subclass_)
    subprocess.run(f"sbatch -J dot_{subclass_} \
        -o ./logs/dot_{subclass_}.log \
        -e ./logs/dot_{subclass_}.log \
        --mem=120G \
        --partition=compute_fat \
        --cpus-per-task=12 \
        --time=150:00:00 \
        ./to_dot.sh {subclass_}",shell=True)
    time.sleep(0.3)
