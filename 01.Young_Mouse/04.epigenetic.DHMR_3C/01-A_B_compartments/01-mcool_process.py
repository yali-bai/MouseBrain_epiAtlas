
import subprocess
import sys
import pandas as pd

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

csv_path = "../../../04.data/02.metainfo/01.RNA/01.Young_Mouse/subclass_corresponding_name.csv"
df = pd.read_csv(csv_path, sep=',')

for subclass_liu in df['subclass_liu']:
    print(f"Processing {subclass_liu}...")

    mcool_file = f'{indir}/{subclass_liu}.Q.100K.mcool'
    result = subprocess.run(['cooler', 'ls', mcool_file], capture_output=True, text=True)
    print(result.stdout)

    c = cooler.Cooler(f"{mcool_file}::/resolutions/100000")
    bins = c.bins()[:]
    print(bins[:5])

    mat_before_balance = c.matrix(balance=False)[:]
    print(mat_before_balance[:5, :5])

    balance_mcool = subprocess.run(
        ['cooler', 'balance', '--force', f'{mcool_file}::/resolutions/100000'],
        stdout=sys.stdout,
        stderr=sys.stderr,
        text=True
    )

    if balance_mcool.returncode != 0:
        print(f"Error: {balance_mcool.stderr}")
    else:
        print("Balance completed successfully")

    bins = c.bins()[:]
    print(bins[:5])

    mat_after_balance = c.matrix(balance=True)[:]
    print(mat_after_balance[:5, :5])

    nan_count = bins['weight'].isna().sum()
    print(f"NAN count: {nan_count}")

    mcool2h5 = subprocess.run([
        'hicConvertFormat', 
        '--matrices', f'{mcool_file}::/resolutions/100000', 
        '--inputFormat', 'cool', 
        '--outputFormat', 'h5', 
        '--outFileName', f'{outdir}/{subclass_liu}_100K.h5'
    ], capture_output=True, text=True)

    print("stdout:", mcool2h5.stdout)
    print("stderr:", mcool2h5.stderr)

    h5_hicPCA = subprocess.run([
        'hicPCA',
        '--matrix', f'{outdir}/{subclass_liu}_100K.h5',
        '--outputFileName', 
        f'{outdir}/{subclass_liu}_100K_PCA_PC1.bedgraph',
        f'{outdir}/{subclass_liu}_100K_PCA_PC2.bedgraph',
        f'{outdir}/{subclass_liu}_100K_PCA_PC3.bedgraph',
        '--format', 'bedgraph',
        '--whichEigenvectors', '1', '2', '3'
    ], capture_output=True, text=True)

    print("stdout:", h5_hicPCA.stdout)
    print("stderr:", h5_hicPCA.stderr)

    pca1_results = pd.read_csv(f'{outdir}/{subclass_liu}_100K_PCA_PC1.bedgraph', sep='\t', header=None)
    pca1_results.columns = ['chromosome', 'start', 'end', 'PC1']

    pca1_filtered = pca1_results[pca1_results['PC1'] != 0].copy()
    pca1_filtered.loc[:, 'compartment'] = ['A' if x > 0 else 'B' for x in pca1_filtered['PC1']]

    compartment_counts_filtered = pca1_filtered['compartment'].value_counts()
    print(compartment_counts_filtered)

    output_path = f'{outdir}/{subclass_liu}_100K_compartments.tsv'
    pca1_filtered.to_csv(output_path, sep='\t', index=False)

    print(f"Completed processing for {subclass_liu}\n")
