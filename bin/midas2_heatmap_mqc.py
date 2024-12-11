#!/usr/bin/env python

import argparse
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        help="Input YAML file containing coverage data"
    )
    parser.add_argument(
        "-s",
        "--samples",
        required=True,
        help="Comma-separated list of all sample names"
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.25,
        help="Minimum fraction of samples a lineage must be present in"
    )
    return parser.parse_args(args)

def main(args):
    # Get all sample names from command line argument
    all_samples = args.samples.split(',')
    
    # Load YAML data
    with open(args.input) as f:
        yaml_data = yaml.safe_load(f)
    
    # Convert to DataFrame
    records = [v for v in yaml_data['data'].values()]
    df = pd.DataFrame(records)
    
    # Identify missing samples
    detected_samples = set(df['sample_name'].unique())
    missing_samples = set(all_samples) - detected_samples
    print(f"Samples with no detectable stats: {', '.join(missing_samples)}")
    
    # Count lineage presence and filter common ones
    lineage_counts = df.groupby('Lineage')['sample_name'].nunique()
    sample_threshold = len(all_samples) * args.threshold
    common_lineages = lineage_counts[lineage_counts > sample_threshold].index
    
    # Filter and pivot data
    filtered_df = df[df['Lineage'].isin(common_lineages)]
    pivot_data = filtered_df.pivot_table(
        index='Lineage',
        columns='Sample',
        values='fraction_covered',
        aggfunc='mean'
    ).reindex(lineage_counts[lineage_counts > sample_threshold].index)
    
    # Add missing samples as columns with zeros
    for sample in missing_samples:
        pivot_data[sample] = 0
    
    # Ensure all samples are included and ordered
    pivot_data = pivot_data.reindex(columns=all_samples, fill_value=0)
    
    # Save data for MultiQC
    pivot_data.round(6).to_csv('midas2_heatmap_mqc.txt', sep='\t')
    
    # Create visualization
    plt.figure(figsize=(15, 10))
    sns.clustermap(
        pivot_data,
        row_cluster=False,
        yticklabels=True if len(pivot_data) <= 30 else False,
        cmap='vlag',
        center=pivot_data.mean().mean(),
        figsize=(16, 18)
    )
    plt.savefig('midas2_heatmap.png')
    plt.close()

if __name__ == "__main__":
    args = parse_args()
    main(args)