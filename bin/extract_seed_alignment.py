#!/usr/bin/env python3
"""
Extract Seed Alignment: Extract and format a sequence's alignment from Rfam seed.

This script reads the Rfam seed alignment file (Rfam.seed or Rfam.seed.gz), finds
the specified Rfam family and sequence, and outputs a Stockholm-format alignment
that can be processed by consensus_projection.py.

This replicates the bpRNA approach of using pre-computed Rfam seed alignments
rather than computing new alignments with cmalign.

Usage:
    extract_seed_alignment.py --seed Rfam.seed.gz --fasta query.fa \\
        --rfam-id RF00001 --seq-id "my_sequence" --output output.sto
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


def open_file(filepath: Path, mode: str = 'rt'):
    """Open a file, handling gzip compression if needed."""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, mode, encoding='utf-8', errors='replace')
    return open(filepath, mode, encoding='utf-8', errors='replace')


def read_fasta_sequence(fasta_path: Path) -> tuple[str, str]:
    """Read a single sequence from a FASTA file, return (id, sequence)."""
    seq_id = None
    sequence = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    break  # Only read first sequence
                seq_id = line[1:].split()[0]
            elif seq_id is not None:
                sequence.append(line.upper())
    
    return seq_id, ''.join(sequence)


def normalize_sequence(seq: str) -> str:
    """Normalize sequence: uppercase, remove gaps, convert T to U."""
    return ''.join(c for c in seq.upper().replace('T', 'U') if c in 'ACGURYMKSWBDHVN')


def extract_family_from_seed(seed_path: Path, rfam_id: str) -> dict:
    """
    Extract a single family's alignment from the Rfam seed file.
    
    Returns a dict with:
        - 'sequences': dict of {seq_name: aligned_seq}
        - 'ss_cons': consensus secondary structure
        - 'rf': reference annotation (if present)
    """
    result = {
        'sequences': {},
        'ss_cons': '',
        'rf': ''
    }
    
    in_target_family = False
    current_accession = None
    
    # For multi-block Stockholm, we need to track parts
    seq_parts = {}
    ss_parts = []
    rf_parts = []
    
    with open_file(seed_path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            
            # Check for new alignment block
            if line.startswith('# STOCKHOLM'):
                in_target_family = False
                current_accession = None
                seq_parts = {}
                ss_parts = []
                rf_parts = []
                continue
            
            # Check for accession
            if line.startswith('#=GF AC'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    # Accession may have version suffix like RF00001.1
                    acc = parts[2].split('.')[0]
                    if acc == rfam_id:
                        in_target_family = True
                        current_accession = acc
                continue
            
            # End of alignment block
            if line.startswith('//'):
                if in_target_family:
                    # Compile the results
                    for seq_name, parts_list in seq_parts.items():
                        result['sequences'][seq_name] = ''.join(parts_list)
                    result['ss_cons'] = ''.join(ss_parts)
                    result['rf'] = ''.join(rf_parts)
                    return result
                continue
            
            # Only process lines if we're in the target family
            if not in_target_family:
                continue
            
            # Parse GC annotations
            if line.startswith('#=GC SS_cons'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    ss_parts.append(parts[2])
                continue
            
            if line.startswith('#=GC RF'):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    rf_parts.append(parts[2])
                continue
            
            # Skip other annotation lines
            if line.startswith('#'):
                continue
            
            # Skip empty lines
            if not line.strip():
                continue
            
            # Parse sequence line
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_name, aligned_seq = parts
                if seq_name not in seq_parts:
                    seq_parts[seq_name] = []
                seq_parts[seq_name].append(aligned_seq)
    
    return result


def find_matching_sequence(family_data: dict, query_seq: str, query_id: str) -> tuple[str, str] | None:
    """
    Find the sequence in the seed alignment that matches the query.
    
    First tries exact match by normalized sequence, then by ID patterns.
    Returns (seq_name, aligned_seq) or None if not found.
    """
    query_normalized = normalize_sequence(query_seq)
    
    # Try exact sequence match
    for seq_name, aligned_seq in family_data['sequences'].items():
        seed_normalized = normalize_sequence(aligned_seq)
        if seed_normalized == query_normalized:
            return seq_name, aligned_seq
    
    # Try partial ID match (query_id might be derived from seed name)
    for seq_name, aligned_seq in family_data['sequences'].items():
        # Common patterns: seq_id might contain accession or partial name
        if query_id in seq_name or seq_name in query_id:
            return seq_name, aligned_seq
    
    # Try subsequence match (query might be a subsequence of seed or vice versa)
    for seq_name, aligned_seq in family_data['sequences'].items():
        seed_normalized = normalize_sequence(aligned_seq)
        if query_normalized in seed_normalized or seed_normalized in query_normalized:
            return seq_name, aligned_seq
    
    return None


def write_stockholm(output_path: Path, seq_name: str, aligned_seq: str, 
                    ss_cons: str, rf_line: str = None):
    """Write a Stockholm-format alignment file."""
    with open(output_path, 'w') as f:
        f.write("# STOCKHOLM 1.0\n\n")
        
        # Write sequence
        f.write(f"{seq_name}  {aligned_seq}\n")
        
        # Write annotations
        if rf_line:
            f.write(f"#=GC RF       {rf_line}\n")
        f.write(f"#=GC SS_cons  {ss_cons}\n")
        
        f.write("//\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequence alignment from Rfam seed file"
    )
    parser.add_argument(
        '--seed', required=True, type=Path,
        help='Path to Rfam.seed or Rfam.seed.gz'
    )
    parser.add_argument(
        '--fasta', required=True, type=Path,
        help='Path to query FASTA file'
    )
    parser.add_argument(
        '--rfam-id', required=True,
        help='Rfam family ID (e.g., RF00001)'
    )
    parser.add_argument(
        '--seq-id', required=True,
        help='Sequence ID for output'
    )
    parser.add_argument(
        '--output', required=True, type=Path,
        help='Output Stockholm file'
    )
    
    args = parser.parse_args()
    
    # Read query sequence
    fasta_id, query_seq = read_fasta_sequence(args.fasta)
    if not query_seq:
        print(f"ERROR: Could not read sequence from {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Query sequence: {args.seq_id} ({len(query_seq)} nt)", file=sys.stderr)
    
    # Extract family from seed
    print(f"Extracting family {args.rfam_id} from seed...", file=sys.stderr)
    family_data = extract_family_from_seed(args.seed, args.rfam_id)
    
    if not family_data['sequences']:
        print(f"ERROR: Family {args.rfam_id} not found in seed file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(family_data['sequences'])} sequences in family", file=sys.stderr)
    
    # Find matching sequence
    match = find_matching_sequence(family_data, query_seq, args.seq_id)
    
    if not match:
        print(f"ERROR: Query sequence not found in seed alignment for {args.rfam_id}", 
              file=sys.stderr)
        print(f"This sequence may not be in the original Rfam seed alignment.", 
              file=sys.stderr)
        print(f"Consider using --alignment_method cmalign instead.", file=sys.stderr)
        sys.exit(1)
    
    seed_name, aligned_seq = match
    print(f"Matched to seed sequence: {seed_name}", file=sys.stderr)
    
    # Write output Stockholm
    write_stockholm(
        args.output,
        args.seq_id,
        aligned_seq,
        family_data['ss_cons'],
        family_data.get('rf')
    )
    
    print(f"Wrote alignment to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
