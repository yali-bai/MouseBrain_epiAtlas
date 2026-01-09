import sys
import subprocess
import datetime
import cooler


subclass_ = sys.argv[1]

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir = ""
# outdir=""

resolution_=10_000

mcool_file=f'{indir}/{subclass_}.Q.10K.mcool'
print(f"Processing {mcool_file}...")

result = subprocess.run(['cooler', 'ls', mcool_file], capture_output=True, text=True)
print(result.stdout)

mcool_file_resolution = f"{mcool_file}::/resolutions/{resolution_}"
print(f'\nUsing {mcool_file_resolution}')

c = cooler.Cooler(mcool_file_resolution)
print(c.bins()[2000:2005])

print(c.matrix(balance=False)[2000:2005, 2000:2005])


balance_mcool = subprocess.run(
    ['cooler', 'balance', '--force', mcool_file_resolution],
    capture_output=True,
    text=True
)

print('balance_mcool finished with returncode:', balance_mcool.returncode)
print('stdout:', balance_mcool.stdout, 'stderr', balance_mcool.stderr, sep='\n')

print(c.bins()[2000:2005])
print(c.matrix(balance=True)[2000:2005, 2000:2005])

print(f"Completed balancing for {subclass_}\n\n")


command=["chromosight", "detect", "--threads", "12", "--min-dist", "20000", "--max-dist", "5000000", mcool_file_resolution, "./output/"+subclass_]


print(subclass_+"\tsubmitted.. with the following command:\n"+" ".join(command))

chromosight_detect = subprocess.run(command, capture_output=True, text=True)

print('chromosight_detect finished with returncode:', chromosight_detect.returncode)
print('stdout:', chromosight_detect.stdout, 'stderr', chromosight_detect.stderr, sep='\n')

print('Finishing now.')