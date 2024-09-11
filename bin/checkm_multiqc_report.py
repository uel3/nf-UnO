#!/usr/bin/env python3

import pandas as pd
import argparse
import yaml

def make_report_yaml(output_file, data_df):
    # Create headers dictionary based on your CheckM columns
    headers = {
        'Bin Id': {
            'title': 'Bin ID'},
        'Marker lineage': {
            'title': 'Marker lineage'},
        '# genomes': {
            'title': '# genomes',
            'format': '{:,.2f}'},
        'Completeness': {
            'title': 'Completeness (%)',
             'description':'Estimated completeness of genome bin based on presence/absence of single copy genes.',
            'format': '{:,.2f}'},
        'Contamination': {
            'title': 'Contamination (%)',
            'description':'Estimated contamination of genome bin based on presence of multiple copies of single copy genes.',
            'format': '{:,.2f}'},
        'Strain heterogeneity': {
            'title':'Strain heterogeneity (%)'},
            'description':'Measure of contamination that can attributed to closely related organisms.',
        'Genome size (bp)': {
            'title': 'Genome size (bp)',
            'description':'Number of nucleotides (including Ns) in the genome bin.',
            'format': '{:,d}'},
        'GC': {
            'title': 'GC (%)',
            'description':'number of G/C nucleotides relative to all nucleotides in the genome bin.',
            'format': '{:,.1f}'},
        '# predicted genes': {
            'title': '# Predicted genes',
            'description':'number of predicted coding sequences (CDS) within the genome bin.',
            'format': '{:,d}'},
        '0':{
            'title':'0',
            'format': '{:,.2f}'},
        '1':{
            'title':'1',
            'format': '{:,.2f}'},
        '2':{
            'title':'2',
            'format': '{:,.2f}'},
        '3':{
            'title':'3',
            'format': '{:,.2f}'},
        '4':{
            'title':'4',
            'format': '{:,.2f}'},
        '5+':{
            'title':'5+',
            'format': '{:,.2f}'}
        # Add more headers as needed
    }

    # Convert the DataFrame to the required format
    data_yaml = data_df.to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': 'checkm_stats',
        'section_name': 'CheckM:Bin Quality Statistics',
        'description': 'Quality assessment of genome bins using CheckM',
        'plot_type': 'table',
        'pconfig': {
            'id': 'checkm_stats',
            'sort_rows': False,
            "scale": False,
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)

def parse_argument():
    parser = argparse.ArgumentParser(prog='create_checkm_report.py')
    parser.add_argument('-i', '--input', metavar='', required=True, help='Specify input CheckM TSV file')
    parser.add_argument('-y', '--yaml', metavar='', required=True, help='Specify output mqc report file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_argument()
    
    # Read the CheckM TSV file
    checkm_df = pd.read_csv(args.input, sep='\t', index_col='Bin Id')
    
    # Select and rename columns if needed
    columns_to_keep = ['Completeness', 'Contamination', 'Strain heterogeneity', 'Genome size (bp)', 'GC', '# predicted genes']
    checkm_df = checkm_df[columns_to_keep]
    
    # Generate report YAML file for MultiQC report
    make_report_yaml(args.yaml, checkm_df)
