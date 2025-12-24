#!/usr/bin/env python3
"""
Parse cmscan tblout output and select the best Rfam hit.

Filters hits based on score and E-value thresholds, selects the best hit
(highest bit score), and outputs result as JSON.
"""

import argparse
import json
import sys
from pathlib import Path


def parse_tblout(tblout_path: Path) -> list[dict]:
    """
    Parse Infernal cmscan tblout format.

    Format (space-separated):
    target_name, target_accession, query_name, query_accession, mdl, mdl_from,
    mdl_to, seq_from, seq_to, strand, trunc, pass, gc, bias, score, E-value,
    inc, description
    """
    hits = []

    with open(tblout_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue

            # Split on whitespace, but description may contain spaces
            parts = line.split()
            if len(parts) < 17:
                continue

            try:
                hit = {
                    'target_name': parts[0],        # Rfam family name
                    'target_accession': parts[1],   # Rfam accession (e.g., RF00001)
                    'query_name': parts[2],         # Sequence name
                    'mdl': parts[4],                # Model type (cm)
                    'mdl_from': int(parts[5]),      # Model start
                    'mdl_to': int(parts[6]),        # Model end
                    'seq_from': int(parts[7]),      # Sequence start
                    'seq_to': int(parts[8]),        # Sequence end
                    'strand': parts[9],             # + or -
                    'trunc': parts[10],             # Truncation info
                    'score': float(parts[14]),      # Bit score
                    'evalue': float(parts[15]),     # E-value
                    'inc': parts[16],               # Inclusion threshold (! or ?)
                }
                hits.append(hit)
            except (ValueError, IndexError) as e:
                sys.stderr.write(f"Warning: Could not parse line: {line}\n")
                continue

    return hits


def filter_hits(hits: list[dict], min_score: float, min_evalue: float) -> list[dict]:
    """Filter hits based on score and E-value thresholds."""
    filtered = []
    for hit in hits:
        if hit['score'] >= min_score and hit['evalue'] <= min_evalue:
            filtered.append(hit)
    return filtered


def select_best_hit(hits: list[dict]) -> dict | None:
    """Select the hit with the highest bit score."""
    if not hits:
        return None
    return max(hits, key=lambda x: x['score'])


def main():
    parser = argparse.ArgumentParser(description='Parse cmscan tblout and select best hit')
    parser.add_argument('--tblout', type=Path, required=True, help='Path to cmscan tblout file')
    parser.add_argument('--seq-id', type=str, required=True, help='Sequence ID')
    parser.add_argument('--min-score', type=float, default=20.0, help='Minimum bit score threshold')
    parser.add_argument('--min-evalue', type=float, default=1e-5, help='Maximum E-value threshold')
    parser.add_argument('--output', type=Path, required=True, help='Output JSON file')

    args = parser.parse_args()

    # Parse tblout file
    hits = parse_tblout(args.tblout)

    # Filter hits
    filtered_hits = filter_hits(hits, args.min_score, args.min_evalue)

    # Select best hit
    best_hit = select_best_hit(filtered_hits)

    # Build result
    if best_hit:
        result = {
            'seq_id': args.seq_id,
            'rfam_id': best_hit['target_accession'],
            'rfam_name': best_hit['target_name'],
            'score': best_hit['score'],
            'evalue': best_hit['evalue'],
            'seq_from': best_hit['seq_from'],
            'seq_to': best_hit['seq_to'],
            'strand': best_hit['strand'],
            'status': 'PASS',
            'total_hits': len(hits),
            'filtered_hits': len(filtered_hits),
        }
    else:
        result = {
            'seq_id': args.seq_id,
            'rfam_id': None,
            'rfam_name': None,
            'score': None,
            'evalue': None,
            'seq_from': None,
            'seq_to': None,
            'strand': None,
            'status': 'NO_HIT',
            'total_hits': len(hits),
            'filtered_hits': 0,
        }

    # Write output
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)

    # Exit with status based on whether we found a hit
    # (used by Nextflow to determine if output is emitted)
    if result['status'] == 'NO_HIT':
        sys.stderr.write(f"No significant hit found for {args.seq_id}\n")


if __name__ == '__main__':
    main()
