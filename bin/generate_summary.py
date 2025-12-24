#!/usr/bin/env python3
"""
Generate summary TSV from projection reports.
"""

import argparse
import json
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Generate summary from projection reports')
    parser.add_argument('--reports', nargs='+', type=Path, help='Projection report JSON files')
    parser.add_argument('--excluded', type=Path, help='Directory with excluded sequence reports')
    parser.add_argument('--output', type=Path, required=True, help='Output TSV file')

    args = parser.parse_args()

    # Collect all reports
    all_reports = []

    # Process successful projections
    if args.reports:
        for report_path in args.reports:
            if report_path.exists() and report_path.suffix == '.json':
                try:
                    with open(report_path, 'r') as f:
                        report = json.load(f)
                        all_reports.append(report)
                except (json.JSONDecodeError, IOError) as e:
                    sys.stderr.write(f"Warning: Could not read {report_path}: {e}\n")

    # Process excluded sequences
    if args.excluded and args.excluded.exists():
        for report_path in args.excluded.glob('*.json'):
            try:
                with open(report_path, 'r') as f:
                    report = json.load(f)
                    all_reports.append(report)
            except (json.JSONDecodeError, IOError) as e:
                sys.stderr.write(f"Warning: Could not read {report_path}: {e}\n")

    # Write TSV
    columns = [
        'seq_id',
        'rfam_id',
        'status',
        'sequence_length',
        'num_base_pairs',
        'num_unpaired',
        'num_insertions',
        'num_noncanonical',
    ]

    with open(args.output, 'w') as f:
        # Header
        f.write('\t'.join(columns) + '\n')

        # Data rows
        for report in all_reports:
            row = [
                report.get('seq_id', ''),
                report.get('rfam_id', '') or '',
                report.get('status', ''),
                str(report.get('sequence_length', '')),
                str(report.get('num_base_pairs', '')),
                str(report.get('num_unpaired', '')),
                str(report.get('num_insertions', '')),
                str(len(report.get('noncanonical_pairs', []))),
            ]
            f.write('\t'.join(row) + '\n')

    sys.stderr.write(f"Summary written: {len(all_reports)} sequences\n")


if __name__ == '__main__':
    main()
