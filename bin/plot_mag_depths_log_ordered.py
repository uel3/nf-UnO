#!/usr/bin/env python
# Originally written by Sabrina Krakau and released under the MIT license.
# See git repository (https://github.com/nf-core/mag) for full license text.
import sys
import argparse
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--bin_depths",
        required=True,
        metavar="FILE",
        help="Bin depths file in TSV format (for one assembly and binning method): bin, sample1_depth, sample2_depth, ....",
    )
    parser.add_argument(
        "-g",
        "--groups",
        required=True,
        metavar="FILE",
        help="File in TSV format containing group information for samples: sample, group",
    )
    parser.add_argument(
        "-o", "--out", required=True, metavar="FILE", type=str, help="Output file."
    )
    return parser.parse_args(args)
def main(args):
    # Load data
    df = pd.read_csv(args.bin_depths, sep="\t", index_col=0)
    groups = pd.read_csv(args.groups, sep="\t", index_col=0, names=["sample", "group"])
    # Count non-zero abundances for each bin
    non_zero_counts = (df > 0).sum(axis=1)
        
    # Sort bins by non-zero abundance count in descending order
    df_sorted = df.loc[non_zero_counts.sort_values(ascending=False).index]
    
    # Log transform the data, adding a small value to handle zeros
    small_value = 1e-6  # can be changed for data scale
    df_log = np.log10(df_sorted + small_value)
    # Prepare colors for group information
    color_map = dict(zip(groups["group"].unique(), sns.color_palette(n_colors=len(groups["group"].unique()))))
   
    # Plot heatmap
    plt.figure(figsize=(12, 10))
    bin_labels = True if len(df_log) <= 30 else False
    g = sns.clustermap(
        df_log,
        row_cluster=False,  # Disable row clustering to maintain our custom order
        yticklabels=bin_labels,
        cmap="vlag",
        center=0,
        col_colors=groups.group.map(color_map),
        figsize=(16, 18)
    )
    g.ax_heatmap.set_xlabel("Samples")
    g.ax_heatmap.set_ylabel("MAGs")
    plt.savefig(args.out)
    plt.close()

    # Prepare data for txt output
    df_output = df_log.round(6)
    # Create the content for the txt file
    txt_content = "Bins\t" + "\t".join(df_output.columns) + "\n"
    for idx, row in df_output.iterrows():
        txt_content += f"{idx}\t" + "\t".join(f"{val}" for val in row) + "\n"
    
    # Save to txt file
    txt_filename = os.path.splitext(args.out)[0] + "_data.txt"
    with open(txt_filename, 'w') as f:
        f.write(txt_content)

if __name__ == "__main__":
    args = parse_args()
    main(args)