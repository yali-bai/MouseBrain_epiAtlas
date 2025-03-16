import sys
import subprocess
import datetime

subclass_ = sys.argv[1]

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

resolution_=10_000

mcool_file=f'{indir}/{subclass_}.Q.10K.mcool::/resolutions/{resolution_}'

command=["chromosight", "detect", "--threads", "12", "--min-dist", "20000", "--max-dist", "5000000", mcool_file, "./output/"+subclass_]


print(subclass_+"\tsubmitted.. with the following command:\n"+" ".join(command))

# Run the command using Popen to capture output in real-time
with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
    for line in proc.stdout:
        with open(f"./logs/{subclass_}.log", "at") as f:
            f.write(str(datetime.datetime.now())+"\n"+line)
    for line in proc.stderr:
        with open(f"./logs/{subclass_}.log", "at") as f:
            f.write(str(datetime.datetime.now())+"\n"+line)

    # Wait for the process to complete
    proc.wait()

# Optionally check if the process was successful
if proc.returncode == 0:
    print(f"{subclass_}:\tCommand executed successfully")
else:
    print(f"{subclass_}:\tError occurred during execution")

