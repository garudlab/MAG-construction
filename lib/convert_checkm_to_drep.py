#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
import re

def reformat_checkm_to_drep(checkm_file, output_file):
    try:
        with open(checkm_file, 'r') as f:
            lines = f.readlines()
        
        header_found = False
        data_lines = []
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('-'):
                continue
            if 'Bin Id' in line and 'Completeness' in line and 'Contamination' in line:
                header_found = True
                continue
            if header_found and line:
                data_lines.append(line)
        
        if not data_lines:
            print("Error: No data found in CheckM output")
            sys.exit(1)
        
        genomes = []
        completeness = []
        contamination = []
        
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 13:
                bin_id = parts[0] + ".fa"
                comp = float(parts[12])
                cont = float(parts[13])
                
                genomes.append(bin_id)
                completeness.append(comp)
                contamination.append(cont)
        
        drep_df = pd.DataFrame({
            'genome': genomes,
            'completeness': completeness,
            'contamination': contamination
        })
        
        drep_df.to_csv(output_file, index=False)
        print(f"Successfully converted {len(drep_df)} genomes to dRep format")
        print(f"Output written to: {output_file}")
        print(f"Sample entries:")
        print(drep_df.head())
        
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Convert CheckM output to dRep format')
    parser.add_argument('-i', '--input', help='CheckM results.tsv file')
    parser.add_argument('-o', '--output', default='genomeInfo.csv', 
                       help='Output file name (default: genomeInfo.csv)')
    
    args = parser.parse_args()
    reformat_checkm_to_drep(args.input, args.output)

if __name__ == "__main__":
    main()
