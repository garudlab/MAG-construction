#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

def reformat_checkm2_to_drep(checkm2_file, output_file):
    try:
        df = pd.read_csv(checkm2_file, sep='\t')

        drep_df = pd.DataFrame({
            'genome': df['Name'] + '.fa',
            'completeness': df['Completeness'],
            'contamination': df['Contamination']
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
    parser = argparse.ArgumentParser(description='Convert CheckM2 output to dRep format')
    parser.add_argument('-i', '--input', help='CheckM2 quality_report.tsv file')
    parser.add_argument('-o', '--output', default='genomeInfo.csv',
                       help='Output file name (default: genomeInfo.csv)')

    args = parser.parse_args()
    reformat_checkm2_to_drep(args.input, args.output)

if __name__ == "__main__":
    main()
