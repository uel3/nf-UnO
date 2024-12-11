#!/usr/bin/env python3
# bin/epi_report.py

import pandas as pd
import numpy as np
import argparse
from typing import Union

# Color codes for terminal output
class Colors:
    # Basic colors
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    
    # Text styles
    BOLD = '\033[1m'

    
    # Reset
    RESET = '\033[0m'

def colorize(text: str, color: str) -> str:
    """Helper function to wrap text in color codes"""
    return f"{color}{text}{Colors.RESET}"

def analyze_bins_and_midas(depth_data, checkm_data, gtdb_data, midas_data, min_samples=2, min_completeness=50, max_contamination=10, midas_threshold=0.25):
    """
    Analyze bins and MIDAS2 data together.
    """
    # Read MAG data
    depths_df = pd.read_csv(depth_data, sep='\t', index_col='bin')
    checkm_df = pd.read_csv(checkm_data, sep='\t', index_col='Sample')
    gtdb_df = pd.read_csv(gtdb_data, sep='\t')
    gtdb_df = gtdb_df.set_index('user_genome')
    
    # Read MIDAS2 data
    midas_df = pd.read_csv(midas_data, sep='\t')
    
    # Process MAG data
    depths_df.index = depths_df.index.str.replace('.fa$', '', regex=True)
    sample_counts = (depths_df > 0).sum(axis=1)
    frequent_bins = sample_counts[sample_counts > min_samples].index
    quality_bins = checkm_df[
        (checkm_df['Completeness'] >= min_completeness) & 
        (checkm_df['Contamination'] <= max_contamination)
    ].index
    valid_bins = set(frequent_bins).intersection(set(quality_bins))
    
    # Process MIDAS2 data
    midas_counts = midas_df.groupby('Lineage')['sample_name'].nunique()
    sample_threshold = len(midas_df['sample_name'].unique()) * midas_threshold
    common_lineages = midas_counts[midas_counts > sample_threshold]
    
    def get_genus_species(taxonomy):
        if pd.isna(taxonomy) or taxonomy == "NA":
            return "NA", "NA"
        taxa = taxonomy.split(';')
        genus = [x for x in taxa if 'g__' in x][0].replace('g__', '')
        species = [x for x in taxa if 's__' in x][0].replace('s__', '')
        return genus, species
    
    summaries = []
    for bin_id in valid_bins:
        presence = depths_df.loc[bin_id]
        present_samples = presence[presence > 0]
        
        taxonomy = gtdb_df.loc[bin_id, 'fastani_taxonomy'] if bin_id in gtdb_df.index else 'NA'
        genus, species = get_genus_species(taxonomy)
        
        midas_match = midas_df[midas_df['Lineage'].str.contains(species, case=False, na=False)] if species != "NA" else pd.DataFrame()
        midas_samples = len(midas_match['sample_name'].unique()) if not midas_match.empty else 0
        midas_coverage = midas_match['fraction_covered'].mean() if not midas_match.empty else 0
        
        summary = {
            'Bin': bin_id,
            'Samples': len(present_samples),
            'Present_in': ', '.join(present_samples.index.tolist()),
            'Completeness': round(checkm_df.loc[bin_id, 'Completeness'], 1),
            'Contamination': round(checkm_df.loc[bin_id, 'Contamination'], 1),
            'Strain_heterogeneity': round(checkm_df.loc[bin_id, 'Strain heterogeneity'], 1),
            'Genome_size': round(checkm_df.loc[bin_id, 'Genome size (bp)']/1000000, 2),
            'GC%': round(checkm_df.loc[bin_id, 'GC'], 1),
            'Taxonomy': f"{genus} {species}",
            'ANI': round(gtdb_df.loc[bin_id, 'fastani_ani'], 2) if bin_id in gtdb_df.index else 'NA',
            'AF': round(gtdb_df.loc[bin_id, 'fastani_af'], 3) if bin_id in gtdb_df.index else 'NA',
            'MIDAS2_samples': midas_samples,
            'MIDAS2_coverage': round(midas_coverage, 3) if midas_coverage > 0 else 'NA'
        }
        summaries.append(summary)
        
    result_df = pd.DataFrame(summaries)
    return result_df.sort_values('Samples', ascending=False).reset_index(drop=True)

def generate_narrative(df, min_completeness=80, max_contamination=10, min_ani=97, min_af=0.8):
    """
    Generate colored narrative with reordered information and quality labels.
    """
    analyzed_bins = df[
        (df['ANI'] != 'NA') &
        (pd.to_numeric(df['ANI'], errors='coerce').notna())
    ]
    
    analyzed_bins = analyzed_bins.sort_values('Samples', ascending=False)
    
    # Calculate summary statistics first
    good_quality_count = 0
    fair_quality_count = 0
    good_quality_midas = 0
    fair_quality_midas = 0
    
    # Pre-calculate quality metrics for summary
    for _, row in analyzed_bins.iterrows():
        completeness_ok = row['Completeness'] >= min_completeness
        contamination_ok = row['Contamination'] <= max_contamination
        ani_ok = row['ANI'] != 'NA' and row['ANI'] >= min_ani
        af_ok = row['AF'] != 'NA' and row['AF'] >= min_af
        
        is_good_quality = all([completeness_ok, contamination_ok, ani_ok, af_ok])
        if is_good_quality:
            good_quality_count += 1
            if row['MIDAS2_samples'] > 0:
                good_quality_midas += 1
        else:
            fair_quality_count += 1
            if row['MIDAS2_samples'] > 0:
                fair_quality_midas += 1
    
    # Start building the narrative with header and summary first
    narrative = f"""{colorize("Summary of MAGs Analysis:", Colors.BOLD)}
{("(Minimum Completeness: " + str(min_completeness) + "%, Maximum Contamination: " + str(max_contamination) + "%, Minimum ANI: " + str(min_ani) + "%, Minimum AF: " + str(min_af) + ")")}

{colorize("Summary Statistics:", Colors.BOLD)}
- Total MAGs analyzed: {str(len(analyzed_bins))}
- Good quality MAGs: {colorize(str(good_quality_count), Colors.GREEN)}
- Fair quality MAGs: {colorize(str(fair_quality_count), Colors.YELLOW)}
- Good quality MAGs with MIDAS2 support: {colorize(str(good_quality_midas), Colors.GREEN)}
- Fair quality MAGs with MIDAS2 support: {colorize(str(fair_quality_midas), Colors.YELLOW)}

{colorize("Most Widely Distributed MAGs:", Colors.BOLD)}"""
    
    # Add detailed MAG information
    for idx, row in analyzed_bins.iterrows():
        # Check each quality metric
        completeness_ok = row['Completeness'] >= min_completeness
        contamination_ok = row['Contamination'] <= max_contamination
        ani_ok = row['ANI'] != 'NA' and row['ANI'] >= min_ani
        af_ok = row['AF'] != 'NA' and row['AF'] >= min_af
        
        # Create list of failed metrics
        failed_metrics = []
        if not completeness_ok:
            failed_metrics.append(f"Completeness below {min_completeness}%")
        if not contamination_ok:
            failed_metrics.append(f"Contamination above {max_contamination}%")
        if not ani_ok:
            failed_metrics.append(f"ANI below {min_ani}%")
        if not af_ok:
            failed_metrics.append(f"AF below {min_af}")
        
        # Determine quality label and message
        is_good_quality = all([completeness_ok, contamination_ok, ani_ok, af_ok])
        if is_good_quality:
            quality_label = colorize("Good quality", Colors.GREEN)
            quality_msg = ""
        else:
            quality_label = colorize("Fair quality", Colors.YELLOW)
            quality_msg = colorize(f" (Failed metrics: {', '.join(failed_metrics)})", Colors.YELLOW)
        
        midas_support = colorize("Found", Colors.GREEN) if row['MIDAS2_samples'] > 0 else colorize("Not found", Colors.RED)
        
        narrative += f"""
{colorize(row['Bin'] + ":", Colors.BOLD)}
- Taxonomy: {row['Taxonomy']} (ANI: {row['ANI']}%, AF: {row['AF']})
- MIDAS2 support: {midas_support} in raw reads (Found in {row['MIDAS2_samples']} samples, Average coverage: {row['MIDAS2_coverage']})
- MAG detection: Present in {row['Samples']} samples ({row['Present_in']})
- Quality:{quality_label}{quality_msg} (Completeness: {row['Completeness']}%, Contamination: {row['Contamination']}%)
- Characteristics: {row['Genome_size']} Mbp genome, {row['GC%']}% GC content, {row['Strain_heterogeneity']}% strain heterogeneity"""
 
    return narrative
def generate_html_narrative(df, min_completeness=80, max_contamination=10, min_ani=97, min_af=0.8):
    """
    Generate HTML-formatted narrative
    """
    html_start = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body { font-family: Arial, sans-serif; max-width: 1200px; margin: 20px auto; padding: 0 20px; }
            h1, h2 { color: #2c3e50; }
            .mag-entry { border: 1px solid #ddd; margin: 10px 0; padding: 15px; border-radius: 5px; }
            .good-quality { color: #27ae60; }
            .fair-quality { color: #f39c12; }
            .midas-found { color: #27ae60; }
            .midas-not-found { color: #c0392b; }
            .failed-metrics { color: #f39c12; font-style: italic; }
            .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin: 20px 0; }
            .parameters { color: #7f8c8d; margin-bottom: 20px; }
        </style>
    </head>
    <body>
    """

    analyzed_bins = df[
        (df['ANI'] != 'NA') &
        (pd.to_numeric(df['ANI'], errors='coerce').notna())
    ]
    
    analyzed_bins = analyzed_bins.sort_values('Samples', ascending=False)
    
    # Track counts for summary
    good_quality_count = 0
    fair_quality_count = 0
    good_quality_midas = 0
    fair_quality_midas = 0
    
    content = f"""
    <h1>Summary of MAGs Analysis</h1>
    <div class="parameters">
    Minimum Completeness: {min_completeness}%, Maximum Contamination: {max_contamination}%, 
    Minimum ANI: {min_ani}%, Minimum AF: {min_af}
    </div>
    """
    
    for idx, row in analyzed_bins.iterrows():
        # Check quality metrics
        completeness_ok = row['Completeness'] >= min_completeness
        contamination_ok = row['Contamination'] <= max_contamination
        ani_ok = row['ANI'] != 'NA' and row['ANI'] >= min_ani
        af_ok = row['AF'] != 'NA' and row['AF'] >= min_af
        
        failed_metrics = []
        if not completeness_ok:
            failed_metrics.append(f"Completeness below {min_completeness}%")
        if not contamination_ok:
            failed_metrics.append(f"Contamination above {max_contamination}%")
        if not ani_ok:
            failed_metrics.append(f"ANI below {min_ani}%")
        if not af_ok:
            failed_metrics.append(f"AF below {min_af}")
        
        is_good_quality = all([completeness_ok, contamination_ok, ani_ok, af_ok])
        if is_good_quality:
            quality_class = "good-quality"
            quality_label = "Good quality"
            quality_msg = ""
            good_quality_count += 1
            if row['MIDAS2_samples'] > 0:
                good_quality_midas += 1
        else:
            quality_class = "fair-quality"
            quality_label = "Fair quality"
            quality_msg = f'<span class="failed-metrics">(Failed metrics: {", ".join(failed_metrics)})</span>'
            fair_quality_count += 1
            if row['MIDAS2_samples'] > 0:
                fair_quality_midas += 1

        midas_class = "midas-found" if row['MIDAS2_samples'] > 0 else "midas-not-found"
        midas_status = "Found" if row['MIDAS2_samples'] > 0 else "Not found"
        
        content += f"""
        <div class="mag-entry">
            <h2>{row['Bin']}</h2>
            <ul>
                <li><strong>Taxonomy:</strong> {row['Taxonomy']} (ANI: {row['ANI']}%, AF: {row['AF']})</li>
                <li><strong>MIDAS2 support:</strong> <span class="{midas_class}">{midas_status}</span> in raw reads 
                    (Found in {row['MIDAS2_samples']} samples, Average coverage: {row['MIDAS2_coverage']})</li>
                <li><strong>MAG detection:</strong> Present in {row['Samples']} samples ({row['Present_in']})</li>
                <li><strong>Quality:</strong> <span class="{quality_class}">{quality_label}</span> {quality_msg}
                    (Completeness: {row['Completeness']}%, Contamination: {row['Contamination']}%)</li>
                <li><strong>Characteristics:</strong> {row['Genome_size']} Mbp genome, {row['GC%']}% GC content, 
                    {row['Strain_heterogeneity']}% strain heterogeneity</li>
            </ul>
        </div>
        """
    
    # Add summary statistics
    summary = f"""
    <div class="summary">
        <h2>Summary Statistics</h2>
        <ul>
            <li>Total MAGs analyzed: {len(analyzed_bins)}</li>
            <li>Good quality MAGs: <span class="good-quality">{good_quality_count}</span></li>
            <li>Fair quality MAGs: <span class="fair-quality">{fair_quality_count}</span></li>
            <li>Good quality MAGs with MIDAS2 support: <span class="good-quality">{good_quality_midas}</span></li>
            <li>Fair quality MAGs with MIDAS2 support: <span class="fair-quality">{fair_quality_midas}</span></li>
        </ul>
    </div>
    """
    
    html_end = """
    </body>
    </html>
    """
    
    return html_start + summary + content + html_end
def main():
    parser = argparse.ArgumentParser(description='Analyze MAG data and generate report')
    parser.add_argument('--depth_data', required=True, help='Path to bin depths summary TSV')
    parser.add_argument('--checkm_data', required=True, help='Path to CheckM stats file')
    parser.add_argument('--gtdb_data', required=True, help='Path to GTDB-Tk summary file')
    parser.add_argument('--midas_data', required=True, help='Path to MIDAS2 abundance file')
    parser.add_argument('--html_output', help='HTML output file name', default='mag_analysis_report.html')
    parser.add_argument('--output', required=True, help='Output file name')
    parser.add_argument('--min_completeness', type=float, default=80, help='Minimum completeness threshold')
    parser.add_argument('--max_contamination', type=float, default=10, help='Maximum contamination threshold')
    parser.add_argument('--min_ani', type=float, default=97, help='Minimum ANI threshold')
    parser.add_argument('--min_af', type=float, default=0.8, help='Minimum AF threshold')

    args = parser.parse_args()

    # Run analysis
    summary_df = analyze_bins_and_midas(
        args.depth_data,
        args.checkm_data,
        args.gtdb_data,
        args.midas_data
    )

    # Generate narrative
    narrative = generate_narrative(
        summary_df,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        min_ani=args.min_ani,
        min_af=args.min_af
    )
    # Generate HTML narrative
    html_narrative = generate_html_narrative(
        summary_df,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        min_ani=args.min_ani,
        min_af=args.min_af
    )

    # Write output
    with open(args.output, 'w') as f:
        f.write(narrative)
    with open(args.html_output, 'w') as f:
        f.write(html_narrative)

if __name__ == '__main__':
    main()

