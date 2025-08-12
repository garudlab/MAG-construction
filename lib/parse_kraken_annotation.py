import pandas as pd
import pickle
from collections import defaultdict, Counter
import argparse

def parse_kraken_line(line):
    """Parse a single Kraken output line"""
    parts = line.strip().split('\t')
    if len(parts) < 5:
        return None
    
    classified = parts[0] == 'C'
    contig_id = parts[1]
    taxid = parts[2]
    length = int(parts[3])
    
    return {
        'classified': classified,
        'contig_id': contig_id,
        'taxid': taxid,
        'length': length
    }

def get_species_name(taxid, taxonomy_dict):
    """
    Get species name for a given taxid using taxonomy mapping
    """
    if str(taxid) not in taxonomy_dict:
        return f"Unknown_taxid_{taxid}"
    
    current_taxid = str(taxid)
    visited = set()
    
    # Traverse up the taxonomy tree looking for species-level classification
    while current_taxid in taxonomy_dict and current_taxid not in visited:
        visited.add(current_taxid)
        info = taxonomy_dict[current_taxid]
        
        # Check if this looks like a species
        if (info['rank_code'].startswith('S') or 
            'species' in info['name'].lower() or
            info['depth'] >= 6):
            return info['name']
        
        # Move to parent
        if info['parent'] != current_taxid:
            current_taxid = info['parent']
        else:
            break
    
    # If no species found, return the original name
    if str(taxid) in taxonomy_dict:
        return taxonomy_dict[str(taxid)]['name']
    
    return f"Unknown_taxid_{taxid}"

def process_kraken_output(kraken_file, contig_mapping_file, taxonomy_mapping_file, 
                          output_file, min_contigs=3, min_fraction=0.1):
    """
    Process Kraken output to get species-level calls per bin
    
    Parameters:
    - kraken_file: path to Kraken output file
    - contig_mapping_file: pickle file with contig to bin mapping
    - taxonomy_mapping_file: pickle file with taxonomy mapping
    - output_file: path to output CSV file
    - min_contigs: minimum number of contigs supporting a species call
    - min_fraction: minimum fraction of bin length supporting a species call
    """
    
    # Load mappings
    print("Loading contig to bin mapping...")
    with open(contig_mapping_file, 'rb') as f:
        contig_to_bin = pickle.load(f)
    
    print("Loading taxonomy mapping...")
    with open(taxonomy_mapping_file, 'rb') as f:
        taxonomy_dict = pickle.load(f)
    
    print(f"Loaded {len(contig_to_bin)} contig mappings and {len(taxonomy_dict)} taxonomy entries")
    
    bin_data = defaultdict(lambda: {
        'total_length': 0,
        'classified_length': 0,
        'taxid_counts': Counter(),
        'taxid_lengths': Counter(),
        'contig_count': 0
    })
    
    # Parse Kraken output
    print("Processing Kraken output...")
    processed_contigs = 0
    mapped_contigs = 0
    
    with open(kraken_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parsed = parse_kraken_line(line)
            if not parsed:
                continue
            
            processed_contigs += 1
            contig_id = parsed['contig_id']
            
            # Get bin name from mapping
            if contig_id not in contig_to_bin:
                continue  # Skip contigs not in our bins
            
            mapped_contigs += 1
            bin_name = contig_to_bin[contig_id]
            length = parsed['length']
            
            bin_data[bin_name]['total_length'] += length
            bin_data[bin_name]['contig_count'] += 1
            
            if parsed['classified'] and parsed['taxid'] != '0':
                bin_data[bin_name]['classified_length'] += length
                bin_data[bin_name]['taxid_counts'][parsed['taxid']] += 1
                bin_data[bin_name]['taxid_lengths'][parsed['taxid']] += length
    
    print(f"Processed {processed_contigs} contigs, {mapped_contigs} mapped to bins")
    
    # Analyze each bin
    results = []
    
    for bin_name, data in bin_data.items():
        total_length = data['total_length']
        classified_length = data['classified_length']
        classification_rate = classified_length / total_length if total_length > 0 else 0
        
        # Find dominant species
        dominant_taxid = None
        dominant_species = "Unclassified"
        dominant_fraction = 0
        dominant_contig_count = 0
        
        if data['taxid_lengths']:
            # Get most abundant taxid by total length
            dominant_taxid = data['taxid_lengths'].most_common(1)[0][0]
            dominant_length = data['taxid_lengths'][dominant_taxid]
            dominant_fraction = dominant_length / total_length
            dominant_contig_count = data['taxid_counts'][dominant_taxid]
            
            # Convert taxid to species name
            dominant_species = get_species_name(dominant_taxid, taxonomy_dict)
        
        # Apply filters
        confidence = "High"
        if dominant_contig_count < min_contigs or dominant_fraction < min_fraction:
            confidence = "Low"
        
        results.append({
            'bin_name': bin_name,
            'species': dominant_species,
            'taxid': dominant_taxid,
            'total_length': total_length,
            'classified_length': classified_length,
            'classification_rate': round(classification_rate, 3),
            'species_fraction': round(dominant_fraction, 3),
            'supporting_contigs': dominant_contig_count,
            'total_contigs': data['contig_count'],
            'confidence': confidence
        })
    
    # Save results
    df = pd.DataFrame(results)
    df = df.sort_values('species_fraction', ascending=False)
    df.to_csv(output_file, index=False)
    
    print(f"Processed {len(results)} bins")
    print(f"Results saved to {output_file}")
    
    # Summary statistics
    high_conf = df[df['confidence'] == 'High']
    print(f"\nSummary:")
    print(f"High confidence classifications: {len(high_conf)}")
    print(f"Low confidence classifications: {len(df) - len(high_conf)}")
    print(f"Mean classification rate: {df['classification_rate'].mean():.3f}")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Convert Kraken contig-level output to bin-level species calls')
    parser.add_argument('--kraken_file', help='Input Kraken output file')
    parser.add_argument('--contig_mapping', help='Pickle file with contig to bin mapping')
    parser.add_argument('--taxonomy_mapping', help='Pickle file with taxonomy mapping')
    parser.add_argument('--output_file', help='Output CSV file')
    parser.add_argument('--min-contigs', type=int, default=3, 
                       help='Minimum contigs supporting species call (default: 3)')
    parser.add_argument('--min-fraction', type=float, default=0.1,
                       help='Minimum fraction of bin length for species call (default: 0.1)')
    
    args = parser.parse_args()
    
    process_kraken_output(args.kraken_file, args.contig_mapping, args.taxonomy_mapping,
                         args.output_file, args.min_contigs, args.min_fraction)

if __name__ == "__main__":
    main()
