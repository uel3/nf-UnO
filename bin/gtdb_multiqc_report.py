#!/usr/bin/env python3
# Script to generate MultiQC-compatible GTDB YAML report

import pandas as pd
import yaml
import os
import re
import argparse

def make_simplified_report_yaml(output_file, data_df):
    """
    Create a simplified YAML file for MultiQC from GTDB-Tk data,
    showing only specific fields and extracting the lowest 2 taxonomic ranks.
    Ensures all dictionary keys are strings for MultiQC compatibility.
    """
    
    # Function to extract the lowest 2 taxonomic levels
    def extract_lowest_taxonomy(taxonomy_string):
        if pd.isna(taxonomy_string):
            return "Unclassified"
        
        # Split the taxonomy string by semicolons
        levels = taxonomy_string.split(';')
        
        # Extract the last two levels (genus and species if available)
        if len(levels) >= 2:
            return ';'.join(levels[-2:])
        else:
            return taxonomy_string

    # Process the DataFrame to extract lowest taxonomy levels
    if 'classification' in data_df.columns:
        data_df['simplified_taxonomy'] = data_df['classification'].apply(extract_lowest_taxonomy)
    
    #if 'fastani_taxonomy' in data_df.columns:
        #data_df['simplified_fastani_taxonomy'] = data_df['fastani_taxonomy'].apply(extract_lowest_taxonomy)
        
    # Select only the requested columns
    columns_to_keep = [
        'user_genome',
        'simplified_taxonomy',  # Our newly created column
        #'fastani_reference',
        'fastani_reference_radius',
        'simplified_fastani_taxonomy',  # Our newly created column
        'fastani_ani',
        'fastani_af',
        'classification_method',
        'note',
        'msa_percent'
    ]
    
    # Filter columns that actually exist in the DataFrame
    existing_columns = [col for col in columns_to_keep if col in data_df.columns]
    simplified_df = data_df[existing_columns].copy()
    
    # Rename the simplified taxonomy columns for display in MultiQC
    if 'simplified_taxonomy' in simplified_df.columns:
        simplified_df.rename(columns={'simplified_taxonomy': 'classification'}, inplace=True)
    
    #if 'simplified_fastani_taxonomy' in simplified_df.columns:
        #simplified_df.rename(columns={'simplified_fastani_taxonomy': 'fastani_taxonomy'}, inplace=True)
    
    # Create headers dictionary based on selected GTDB-Tk columns
    headers = {
        'user_genome': {
            'title': 'Genome ID',
            'description': 'Genome bin identifier'
        },
        'classification': {
            'title': 'Classification',
            'description': 'Lowest 2 taxonomic levels'
        },
        'fastani_reference': {
            'title': 'FastANI Reference',
            'description': 'Reference genome used for FastANI comparison'
        },
        'fastani_reference_radius': {
            'title': 'FastANI Ref Radius',
            'description': 'Radius used for FastANI reference comparison',
            'format': '{:,.1f}'
        },
        #'fastani_taxonomy': {
            #'title': 'FastANI Taxonomy',
            #'description': 'Lowest 2 taxonomic levels of FastANI reference'
        #},
        'fastani_ani': {
            'title': 'ANI',
            'description': 'Average Nucleotide Identity value from FastANI',
            'format': '{:,.2f}'
        },
        'fastani_af': {
            'title': 'AF',
            'description': 'Alignment fraction from FastANI',
            'format': '{:,.3f}'
        },
        'classification_method': {
            'title': 'Method',
            'description': 'Method used for taxonomic classification'
        },
        'note': {
            'title': 'Note',
            'description': 'Additional information about the classification'
        },
        'msa_percent': {
            'title': 'MSA %',
            'description': 'Percentage of markers in the MSA',
            'format': '{:,.2f}'
        }
    }
    
    # Keep only headers for columns that exist in the simplified DataFrame
    existing_headers = {k: v for k, v in headers.items() if k in simplified_df.columns}
    
    # Convert the DataFrame to the required format with string keys for index values
    data_yaml = {}
    for idx, row in simplified_df.iterrows():
        # Use user_genome as the key if available, otherwise use string index
        if 'user_genome' in row:
            key = str(row['user_genome'])
        else:
            key = f"bin_{str(idx)}"
        
        # Create a copy of the row as a dictionary
        row_dict = row.to_dict()
        
        # Add the data to the YAML dictionary
        data_yaml[key] = row_dict
    
    # Create the full YAML dictionary with names matching your MultiQC config
    yaml_dict = {
        'id': 'gtdb_summary',  # Match the ID in your MultiQC config
        'section_name': 'GTDB Taxonomic Classifications',  # Match the section_name in your config
        'description': 'Taxonomic Classification of genome bins using GTDB',  # Match the description in your config
        'plot_type': 'table',
        'pconfig': {
            'id': 'gtdb_summary',  # Match the ID in your MultiQC config
            'sort_rows': False,
            'scale': False,
        },
        'headers': existing_headers,
        'data': data_yaml
    }
    
    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)
    
    return len(simplified_df)

def parse_args():
    parser = argparse.ArgumentParser(description='Generate MultiQC-compatible YAML from GTDB-Tk summary file')
    parser.add_argument('-i', '--input', required=True, help='Path to GTDB-Tk summary TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output YAML file path (e.g., sample_id_gtdb_mqc.yml)')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    try:
        # Read the GTDB-Tk TSV file
        gtdbtk_df = pd.read_csv(args.input, sep='\t')
        
        # Generate simplified report YAML file for MultiQC report
        make_simplified_report_yaml(args.output, gtdbtk_df)
        
        # Verify the YAML can be read back
        with open(args.output, 'r') as f:
            yaml.safe_load(f)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())