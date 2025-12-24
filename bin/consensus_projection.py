#!/usr/bin/env python3
"""
Consensus Projection: Project Rfam consensus secondary structure onto a query sequence.

This script reads a Stockholm alignment (from cmalign) and projects the consensus
secondary structure (#=GC SS_cons) onto the query sequence, handling gaps and
insertions appropriately.

Projection Rules:
1. Match state (query aligned to consensus): Inherit SS_cons character
2. Insertion state (query has residue, consensus has gap): Assign '.' (unpaired)
3. Deletion state (query has gap, consensus has residue): Skip (no residue to annotate)
4. Broken pairs: If one partner is deleted, convert remaining bracket to '.'
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import NamedTuple


class AlignmentData(NamedTuple):
    """Parsed Stockholm alignment data."""
    seq_id: str
    aligned_seq: str
    ss_cons: str
    rf_line: str | None  # Reference annotation line (if present)


# Canonical Watson-Crick and wobble base pairs
CANONICAL_PAIRS = {
    ('G', 'C'), ('C', 'G'),
    ('A', 'U'), ('U', 'A'),
    ('G', 'U'), ('U', 'G'),
    # Also handle T as U
    ('A', 'T'), ('T', 'A'),
    ('G', 'T'), ('T', 'G'),
}

# Bracket pairs for structure notation
OPENING_BRACKETS = '(<[{'
CLOSING_BRACKETS = ')>]}'
BRACKET_PAIRS = dict(zip(OPENING_BRACKETS, CLOSING_BRACKETS))
REVERSE_BRACKET_PAIRS = dict(zip(CLOSING_BRACKETS, OPENING_BRACKETS))


def parse_stockholm(sto_path: Path) -> AlignmentData:
    """
    Parse a Stockholm alignment file from cmalign.

    Expected format:
    # STOCKHOLM 1.0
    #=GC SS_cons  <<<..>>>
    #=GC RF       xxxx.xxx
    seq_name      ACGU-CGU
    //
    """
    aligned_seq = None
    ss_cons = None
    rf_line = None
    seq_id = None

    # For multi-block Stockholm files, we need to concatenate
    seq_parts = []
    ss_parts = []
    rf_parts = []

    with open(sto_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')

            # Skip comments and empty lines (except GC lines)
            if line.startswith('# STOCKHOLM') or line.startswith('//') or not line:
                continue

            # Parse GC (per-column) annotations
            if line.startswith('#=GC SS_cons'):
                # Format: #=GC SS_cons   <structure>
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

            # Parse sequence line
            # Format: seq_name   ALIGNED_SEQUENCE
            parts = line.split(None, 1)
            if len(parts) == 2:
                if seq_id is None:
                    seq_id = parts[0]
                if parts[0] == seq_id:
                    seq_parts.append(parts[1].replace(' ', ''))

    # Concatenate parts
    aligned_seq = ''.join(seq_parts) if seq_parts else None
    ss_cons = ''.join(ss_parts) if ss_parts else None
    rf_line = ''.join(rf_parts) if rf_parts else None

    if aligned_seq is None or ss_cons is None:
        raise ValueError(f"Could not parse Stockholm file: missing sequence or SS_cons")

    if len(aligned_seq) != len(ss_cons):
        raise ValueError(
            f"Alignment length mismatch: sequence={len(aligned_seq)}, SS_cons={len(ss_cons)}"
        )

    return AlignmentData(
        seq_id=seq_id,
        aligned_seq=aligned_seq,
        ss_cons=ss_cons,
        rf_line=rf_line,
    )


def normalize_ss_cons(ss_cons: str) -> str:
    """
    Normalize SS_cons notation.

    Converts various bracket notations to standard dot-bracket:
    - Replaces '_', '-', ':', ',' with '.' (unpaired)
    - Keeps brackets: (), <>, [], {}
    - Keeps letters for pseudoknots: A-Z (opening), a-z (closing)
    """
    result = []
    for char in ss_cons:
        if char in OPENING_BRACKETS or char in CLOSING_BRACKETS:
            result.append(char)
        elif char == '.':
            result.append('.')
        elif char in '_-:,~':
            # These represent unpaired positions in Rfam notation
            result.append('.')
        elif char.isalpha():
            # Letters represent pseudoknot pairs (uppercase=open, lowercase=close)
            result.append(char)
        else:
            # Unknown character, treat as unpaired
            result.append('.')
    return ''.join(result)


def find_pairs(structure: str) -> dict[int, int]:
    """
    Find base pair positions in a bracket notation structure.

    Handles:
    - Bracket pairs: (), <>, [], {}
    - Letter-based pseudoknots: A-Z (opening) pairs with a-z (closing)

    Returns a dictionary mapping each paired position to its partner.
    """
    pairs = {}
    stacks = {bracket: [] for bracket in OPENING_BRACKETS}
    # Add stacks for letter-based pseudoknots (A-Z)
    letter_stacks = {chr(c): [] for c in range(ord('A'), ord('Z') + 1)}

    for i, char in enumerate(structure):
        if char in OPENING_BRACKETS:
            stacks[char].append(i)
        elif char in CLOSING_BRACKETS:
            opening = REVERSE_BRACKET_PAIRS[char]
            if stacks[opening]:
                j = stacks[opening].pop()
                pairs[i] = j
                pairs[j] = i
        elif char.isupper():
            # Uppercase letter = opening pseudoknot
            letter_stacks[char].append(i)
        elif char.islower():
            # Lowercase letter = closing pseudoknot
            upper = char.upper()
            if letter_stacks[upper]:
                j = letter_stacks[upper].pop()
                pairs[i] = j
                pairs[j] = i

    return pairs


def project_structure(aligned_seq: str, ss_cons: str) -> tuple[str, str, list[tuple[int, int]]]:
    """
    Project consensus structure onto query sequence.

    Args:
        aligned_seq: Aligned query sequence (with gaps)
        ss_cons: Consensus secondary structure

    Returns:
        - Ungapped sequence
        - Projected structure (same length as ungapped sequence)
        - List of (query_pos, aln_col) for insertions
    """
    # Normalize the consensus structure
    ss_cons = normalize_ss_cons(ss_cons)

    # Find pairs in the alignment coordinates
    aln_pairs = find_pairs(ss_cons)

    # Track which alignment columns are match states vs insertions
    # In cmalign output: lowercase = insertion, uppercase = match
    # Gap in RF line (if present) also indicates insertion

    projected_structure = []
    ungapped_sequence = []
    insertions = []  # (query_pos, alignment_col)

    # Map from alignment column to query position
    aln_to_query = {}
    query_pos = 0

    for aln_col, (seq_char, ss_char) in enumerate(zip(aligned_seq, ss_cons)):
        if seq_char not in '.-':
            # Query has a residue at this position
            aln_to_query[aln_col] = query_pos

            # Check if this is an insertion (lowercase in cmalign output)
            is_insertion = seq_char.islower()

            if is_insertion:
                # Insertion: assign unpaired
                projected_structure.append('.')
                insertions.append((query_pos, aln_col))
            else:
                # Match state: inherit structure
                projected_structure.append(ss_char)

            ungapped_sequence.append(seq_char.upper())
            query_pos += 1

    # Now fix broken pairs
    # A pair is broken if one partner was deleted (gap in query)
    projected_list = list(projected_structure)

    for aln_col, ss_char in enumerate(ss_cons):
        # Check both bracket pairs and letter-based pseudoknots
        is_paired_char = (ss_char in OPENING_BRACKETS or 
                         ss_char in CLOSING_BRACKETS or 
                         ss_char.isalpha())
        if is_paired_char:
            # Check if this position exists in query
            if aln_col in aln_to_query:
                query_idx = aln_to_query[aln_col]
                # Check if partner exists
                if aln_col in aln_pairs:
                    partner_aln_col = aln_pairs[aln_col]
                    if partner_aln_col not in aln_to_query:
                        # Partner was deleted - break the pair
                        projected_list[query_idx] = '.'

    return ''.join(ungapped_sequence), ''.join(projected_list), insertions


def validate_base_pairs(
    sequence: str,
    structure: str,
    mode: str = 'keep'
) -> tuple[str, list[dict]]:
    """
    Validate base pairs for canonical pairing.

    Args:
        sequence: RNA sequence
        structure: Dot-bracket structure
        mode: 'keep' (preserve all), 'remove' (break non-canonical), 'flag' (keep but report)

    Returns:
        - Possibly modified structure
        - List of non-canonical pairs
    """
    pairs = find_pairs(structure)
    structure_list = list(structure)
    noncanonical = []

    # Only check each pair once (when we see the opening bracket or uppercase letter)
    checked = set()

    for i, char in enumerate(structure):
        # Check opening brackets and uppercase letters (pseudoknot openers)
        is_opener = char in OPENING_BRACKETS or char.isupper()
        if is_opener and i in pairs:
            j = pairs[i]
            if (i, j) in checked or (j, i) in checked:
                continue
            checked.add((i, j))

            bp = (sequence[i].upper(), sequence[j].upper())
            if bp not in CANONICAL_PAIRS:
                noncanonical.append({
                    'pos_i': i,
                    'pos_j': j,
                    'base_i': sequence[i],
                    'base_j': sequence[j],
                })
                if mode == 'remove':
                    structure_list[i] = '.'
                    structure_list[j] = '.'

    return ''.join(structure_list), noncanonical


def standardize_brackets(structure: str) -> str:
    """
    Convert all bracket types to standard parentheses ().

    Converts: <>, [], {} -> ()
    Preserves: Letter-based pseudoknots (A...a, B...b, etc.) and dots
    """
    result = []
    for char in structure:
        if char in OPENING_BRACKETS:
            result.append('(')
        elif char in CLOSING_BRACKETS:
            result.append(')')
        else:
            # Keep letters (pseudoknots) and dots as-is
            result.append(char)
    return ''.join(result)


def write_dotbracket(
    output_path: Path,
    seq_id: str,
    rfam_id: str,
    sequence: str,
    structure: str,
    score: float | None = None,
    status: str = 'PASS',
    preserve_pseudoknots: bool = False
):
    """Write sequence and structure in dot-bracket format."""
    # Optionally standardize all brackets to ()
    if not preserve_pseudoknots:
        structure = standardize_brackets(structure)

    header_parts = [seq_id, rfam_id]
    if score is not None:
        header_parts.append(f"score={score:.1f}")
    header_parts.append(f"status={status}")

    with open(output_path, 'w') as f:
        f.write(f">{' '.join(header_parts)}\n")
        f.write(f"{sequence}\n")
        f.write(f"{structure}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Project Rfam consensus structure onto query sequence'
    )
    parser.add_argument('--input', type=Path, required=True, help='Stockholm alignment file')
    parser.add_argument('--output', type=Path, required=True, help='Output dot-bracket file')
    parser.add_argument('--report', type=Path, required=True, help='Output JSON report')
    parser.add_argument('--seq-id', type=str, required=True, help='Sequence ID')
    parser.add_argument('--rfam-id', type=str, required=True, help='Rfam family ID')
    parser.add_argument(
        '--noncanonical',
        choices=['keep', 'remove', 'flag'],
        default='keep',
        help='How to handle non-canonical base pairs'
    )
    parser.add_argument(
        '--preserve-pseudoknots',
        action='store_true',
        help='Preserve pseudoknot bracket notation (<>, [], {}) instead of converting to ()'
    )
    parser.add_argument('--score', type=float, help='cmscan bit score (for header)')

    args = parser.parse_args()

    # Parse Stockholm alignment
    try:
        aln_data = parse_stockholm(args.input)
    except Exception as e:
        sys.stderr.write(f"Error parsing Stockholm file: {e}\n")
        sys.exit(1)

    # Project structure
    sequence, structure, insertions = project_structure(
        aln_data.aligned_seq,
        aln_data.ss_cons
    )

    # Validate base pairs
    structure, noncanonical = validate_base_pairs(
        sequence,
        structure,
        mode=args.noncanonical
    )

    # Count structural features
    pairs = find_pairs(structure)
    num_pairs = len(pairs) // 2
    num_unpaired = structure.count('.')

    # Write dot-bracket output
    write_dotbracket(
        args.output,
        args.seq_id,
        args.rfam_id,
        sequence,
        structure,
        score=args.score,
        status='PASS',
        preserve_pseudoknots=args.preserve_pseudoknots
    )

    # Write JSON report
    report = {
        'seq_id': args.seq_id,
        'rfam_id': args.rfam_id,
        'sequence_length': len(sequence),
        'num_base_pairs': num_pairs,
        'num_unpaired': num_unpaired,
        'num_insertions': len(insertions),
        'noncanonical_pairs': noncanonical,
        'noncanonical_mode': args.noncanonical,
        'status': 'PASS',
    }

    with open(args.report, 'w') as f:
        json.dump(report, f, indent=2)

    # Summary to stderr
    sys.stderr.write(
        f"Projected structure for {args.seq_id}: "
        f"{len(sequence)} nt, {num_pairs} pairs, {len(insertions)} insertions\n"
    )


if __name__ == '__main__':
    main()
