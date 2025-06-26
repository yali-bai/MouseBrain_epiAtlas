import subprocess
import pandas as pd

with open("../../../04.data/04.config_files/subclass_order_for_integration_with_zeng.txt", "rt") as f:
    subclass_list=f.read().split("\n")[:-1]
subclass_list=[i for i in subclass_list if "NN" not in i]

csv_path = "../../../04.data/02.metainfo/01.RNA/01.Young_Mouse/subclass_corresponding_name.csv"
df = pd.read_csv(csv_path, sep=',')
df = df.iloc[:len(subclass_list),:]


from concurrent.futures import ThreadPoolExecutor
## for salloc and srun --pty bash!!! change to sbatch if on mgt!!

# Define the function to run the subprocess
def run_subprocess(subclass_):
    subprocess.run(f"bash ./to_dot.sh {subclass_}", shell=True)

# Using ThreadPoolExecutor to run subprocess calls in parallel
with ThreadPoolExecutor() as executor:
    # Submit tasks to the executor
    executor.map(run_subprocess, df["subclass_liu"])