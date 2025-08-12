from pathlib import Path
from Bio import SeqIO
import pickle
import pandas as pd
from collections import defaultdict
import argparse

def build_contig_to_bin_mapping(genomes_dir):
    """
    Build mapping from contig names to bin names by reading all .fa files
    """
    genomes_path = Path(genomes_dir)
    contig_to_bin = {}
    
    print(f"Scanning {genomes_path} for .fa files...")
    
    fa_files = list(genomes_path.glob("*.fa")) + list(genomes_path.glob("*.fasta"))
    
    if not fa_files:
        print("No .fa or .fasta files found!")
        return {}
    
    print(f"Found {len(fa_files)} genome files")
    
    for fa_file in fa_files:
        bin_name = fa_file.stem  # filename without extension
        print(f"Processing {bin_name}...")
        
        try:
            for record in SeqIO.parse(fa_file, "fasta"):
                contig_name = record.id
                contig_to_bin[contig_name] = bin_name
        except Exception as e:
            print(f"Error processing {fa_file}: {e}")
            continue
    
    print(f"Built mapping for {len(contig_to_bin)} contigs across {len(fa_files)} bins")
    return contig_to_bin

def build_taxonomy_mapping(taxonomy_file):
    """
    Build taxonomy mapping from ktaxonomy.tsv
    Format appears to be: taxid | parent_taxid | rank_code | depth | name
    """
    print(f"Building taxonomy mapping from {taxonomy_file}...")
    
    taxid_to_info = {}
    
    with open(taxonomy_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = [p.strip() for p in line.strip().split('|')]
            if len(parts) >= 5:
                taxid = parts[0]
                parent_taxid = parts[1] 
                rank_code = parts[2]
                depth = parts[3]
                name = parts[4]
                
                taxid_to_info[taxid] = {
                    'parent': parent_taxid,
                    'rank_code': rank_code,
                    'depth': int(depth) if depth.isdigit() else 0,
                    'name': name
                }
    
    print(f"Built taxonomy mapping for {len(taxid_to_info)} taxa")
    return taxid_to_info

def get_species_name(taxid, taxonomy_dict):
    """
    Get species name for a given taxid by traversing up the taxonomy
    """
    if str(taxid) not in taxonomy_dict:
        return f"Unknown_taxid_{taxid}"
    
    current_taxid = str(taxid)
    visited = set()
    
    # Traverse up the taxonomy tree looking for species-level classification
    while current_taxid in taxonomy_dict and current_taxid not in visited:
        visited.add(current_taxid)
        info = taxonomy_dict[current_taxid]
        
        # Check if this looks like a species (you may need to adjust these conditions)
        if (info['rank_code'].startswith('S') or 
            'species' in info['name'].lower() or
            info['depth'] >= 6):  # Species typically at depth 6+
            return info['name']
        
        # Move to parent
        if info['parent'] != current_taxid:  # avoid infinite loops
            current_taxid = info['parent']
        else:
            break
    
    # If no species found, return the original name or genus
    if str(taxid) in taxonomy_dict:
        return taxonomy_dict[str(taxid)]['name']
    
    return f"Unknown_taxid_{taxid}"

def main():
    parser = argparse.ArgumentParser(description='Build contig-to-bin mapping and taxonomy lookup')
    parser.add_argument('--genomes_dir', help='Directory containing .fa genome files')
    parser.add_argument('--taxonomy_file', help='ktaxonomy.tsv file', default='/u/project/ngarud/Garud_lab/averster/databases/kraken/ktaxonomy.tsv')
    parser.add_argument('--contig-mapping-out', default='/u/project/ngarud/Garud_lab/averster/mag_results/04_dRep/contig_to_bin.pkl',
                       help='Output pickle file for contig mapping')
    parser.add_argument('--taxonomy-mapping-out', default='/u/project/ngarud/Garud_lab/averster/mag_results/04_dRep/taxonomy_mapping.pkl',
                       help='Output pickle file for taxonomy mapping')
    
    args = parser.parse_args()
    
    # Build contig to bin mapping
    print("Building contig to bin mapping...")
    contig_to_bin = build_contig_to_bin_mapping(args.genomes_dir)
    
    if contig_to_bin:
        with open(args.contig_mapping_out, 'wb') as f:
            pickle.dump(contig_to_bin, f)
        print(f"Contig mapping saved to {args.contig_mapping_out}")
        
        # Show some examples
        print("\nExample mappings:")
        for i, (contig, bin_name) in enumerate(contig_to_bin.items()):
            if i < 5:
                print(f"  {contig} -> {bin_name}")
            else:
                break
    
    # Build taxonomy mapping
    print("\nBuilding taxonomy mapping...")
    taxonomy_dict = build_taxonomy_mapping(args.taxonomy_file)
    
    if taxonomy_dict:
        with open(args.taxonomy_mapping_out, 'wb') as f:
            pickle.dump(taxonomy_dict, f)
        print(f"Taxonomy mapping saved to {args.taxonomy_mapping_out}")
        
        # Test a few species lookups
        print("\nTesting species lookups:")
        test_taxids = ['562', '1239', '2608920']  # Some common bacterial taxids
        for taxid in test_taxids:
            species = get_species_name(taxid, taxonomy_dict)
            print(f"  taxid {taxid} -> {species}")

if __name__ == "__main__":
    main()
