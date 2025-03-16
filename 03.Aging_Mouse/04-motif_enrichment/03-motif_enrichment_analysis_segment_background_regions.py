import pandas as pd
import tempfile
import subprocess
import os
import sys
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from intervaltree import Interval, IntervalTree

def load_background_regions(file_path):
    """
    Load background regions from a file, assuming format: chromosome, start, end
    :param file_path: Path to the background regions file
    :return: DataFrame containing the background regions
    """
    background_df = pd.read_csv(file_path, sep='\t', header=None, names=['chromosome', 'start', 'end'])
    return background_df

def remove_overlapping_regions(background_df, all_subclass_regions):
    """
    Remove background regions that overlap with any subclass foreground regions
    :param background_df: DataFrame of background regions
    :param all_subclass_regions: Merged list of all subclass foreground regions
    :return: DataFrame of background regions with overlaps removed
    """
    # Create an interval tree for subclass regions
    subclass_tree = IntervalTree()
    for chrom, start, end in all_subclass_regions:
        start = int(start)
        end = int(end)
        if start < end and pd.notnull(start) and pd.notnull(end):  # Check if the interval is valid
            subclass_tree.add(Interval(start, end, chrom))
    
    non_overlapping_regions = []

    for _, background_row in background_df.iterrows():
        background_chrom = background_row['chromosome']
        background_start = int(background_row['start'])
        background_end = int(background_row['end'])
        
        # Check if the background region overlaps with any subclass region in the same chromosome
        overlaps = subclass_tree[background_start:background_end]
        overlaps = [overlap for overlap in overlaps if overlap.data == background_chrom]

        if not overlaps:
            non_overlapping_regions.append(background_row)
    
    return pd.DataFrame(non_overlapping_regions, columns=['chromosome', 'start', 'end'])

def annotate_regions(df, motif_df):
    """
    Annotate regions with motif data using bedtools.
    :param df: DataFrame of regions to annotate
    :param motif_df: DataFrame of motifs to use for annotation
    :return: DataFrame with annotated regions
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed')
    df.to_csv(temp_file.name, sep='\t', header=False, index=False)
    
    motif_temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed')
    motif_df.to_csv(motif_temp_file.name, sep='\t', header=False, index=False)
    
    intersect_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed').name
    command = f'bedtools intersect -a {temp_file.name} -b {motif_temp_file.name} -r -wa -wb > {intersect_file}'
    subprocess.run(command, shell=True, check=True)
    
    intersect_data = pd.read_csv(intersect_file, sep='\t', header=None, names=[
        'chromosome', 'start', 'end', 'motif_chromosome', 'motif_start', 'motif_end', 'motif', 'score', 'pvalue', 'strand'
    ])
    
    os.remove(temp_file.name)
    os.remove(motif_temp_file.name)
    os.remove(intersect_file)
    
    return intersect_data.drop_duplicates(subset=['chromosome', 'start', 'end'])

def calculate_enrichment(subclass_df, background_df, subclass_total_regions, background_total_regions):
    """
    Calculate motif enrichment and perform Fisher's Exact test.
    :param subclass_df: Motif data for subclass regions
    :param background_df: Motif data for background regions
    :param subclass_total_regions: Total number of subclass regions
    :param background_total_regions: Total number of background regions
    :return: DataFrame with enrichment results
    """
    subclass_motif_counts = subclass_df['motif'].value_counts()
    background_motif_counts = background_df['motif'].value_counts()
    
    all_motifs_set = set(subclass_motif_counts.index) | set(background_motif_counts.index)
    all_motifs = pd.DataFrame({
        'motif': list(all_motifs_set),
        'count_subclass': [subclass_motif_counts.get(motif, 0) for motif in all_motifs_set],
        'count_background': [background_motif_counts.get(motif, 0) for motif in all_motifs_set]
    })
    
    results = []
    for _, row in all_motifs.iterrows():
        motif = row['motif']
        count_subclass = row['count_subclass']
        count_background = row['count_background']
        if count_subclass == 0 and count_background == 0:
            continue
        oddsratio, p_value = fisher_exact([
            [count_subclass, subclass_total_regions - count_subclass],
            [count_background, background_total_regions - count_background]
        ])
        results.append({
            'motif': motif,
            'count_subclass': count_subclass,
            'total_subclass_minus_count': subclass_total_regions - count_subclass,
            'count_background': count_background,
            'total_background_minus_count': background_total_regions - count_background,
            'oddsratio': oddsratio,
            'p_value': p_value
        })
    
    results_df = pd.DataFrame(results)
    results_df['adjusted_p_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
    return results_df

def main(input_csv, motif_tsv_gz, output_csv, region_type, background_file):
    """
    Main function to perform enrichment analysis.
    :param input_csv: Input file containing subclass region data
    :param motif_tsv_gz: Motif database file
    :param output_csv: Output file path for results
    :param region_type: Type of region to analyze
    :param background_file: Path to the background region file
    """
    df = pd.read_csv(input_csv, index_col=0)
    df = df[['Subclass', region_type]]
    df['parsed_regions'] = df[region_type].apply(lambda x: [region.split('_') for region in eval(x)])
    
    motif_data = pd.read_csv(motif_tsv_gz, sep='\t', compression='gzip')
    motif_data.columns = ['chromosome', 'start', 'end', 'motif', 'score', 'pvalue', 'strand']
    
    background_df = load_background_regions(background_file)  # Load background regions
    background_df = background_df.dropna(subset=['start', 'end'])
    background_df = background_df[background_df['start'] < background_df['end']]
    
    # Merge all subclass regions into one list
    all_subclass_regions = []
    for _, row in df.iterrows():
        all_subclass_regions.extend(row['parsed_regions'])
    
    all_subclass_regions = [
    (chrom, start, end) for chrom, start, end in all_subclass_regions
    if start < end and pd.notnull(start) and pd.notnull(end)
    ]
    
    # Remove overlapping regions with any subclass foreground regions
    background_df = remove_overlapping_regions(background_df, all_subclass_regions)
    
    # Calculate total number of background regions
    background_total_regions = len(background_df)  # Total number of background regions after removing overlaps
    
    # Annotate background regions
    background_df_annotated = annotate_regions(background_df, motif_data)
    
    all_results = []
    
    for subclass in df['Subclass'].unique():
        subclass_row = df[df['Subclass'] == subclass].iloc[0]
        subclass_df = annotate_regions(pd.DataFrame(subclass_row['parsed_regions'], columns=['chromosome', 'start', 'end']), motif_data)
        subclass_total_regions = len(subclass_row['parsed_regions'])
        
        # Perform enrichment analysis
        enrichment_results = calculate_enrichment(subclass_df, background_df_annotated, subclass_total_regions, background_total_regions)
        enrichment_results['subclass'] = subclass
        all_results.append(enrichment_results)
    
    final_results = pd.concat(all_results, ignore_index=True)
    final_results.to_csv(output_csv, index=False)
    print(f"Enrichment analysis results saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python 03-motif_enrichment_analysis_segment_background_regions.py <input_csv> <motif_tsv_gz> <output_csv> <region_type> <background_file>")
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
