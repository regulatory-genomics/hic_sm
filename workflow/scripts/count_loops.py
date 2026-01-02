#!/usr/bin/env python3
"""
Count loops from mustache loop detection output files across multiple resolutions.

Reads loop TSV files (one per resolution), counts the number of loops in each,
and generates a summary statistics file with counts per resolution and total.
"""

import argparse
import sys
import os
from pathlib import Path
import pandas as pd


def count_loops_in_file(loop_file):
    """
    Count the number of loops in a loop TSV file.
    
    Args:
        loop_file: Path to loop TSV file
    
    Returns:
        int: Number of loops (rows excluding header)
    """
    if not os.path.exists(loop_file):
        return 0
    
    try:
        # Read the file and count non-header rows
        df = pd.read_csv(loop_file, sep='\t')
        return len(df)
    except pd.errors.EmptyDataError:
        return 0
    except Exception as e:
        print(f"Warning: Error reading {loop_file}: {e}", file=sys.stderr)
        return 0


def extract_resolution_from_filename(filename):
    """
    Extract resolution from filename pattern: {prefix}_{resolution}_loops.tsv
    
    Args:
        filename: Path or filename string
    
    Returns:
        str: Resolution value or filename if pattern doesn't match
    """
    filename = Path(filename).name
    # Pattern: {library}.{filter_name}.{resolution}_loops.tsv
    if '_loops.tsv' in filename:
        parts = filename.replace('_loops.tsv', '').split('.')
        # Resolution is typically the last part before _loops
        # Try to find numeric resolution
        for part in reversed(parts):
            if part.isdigit():
                return part
        # If no numeric found, return the part before _loops
        return filename.replace('_loops.tsv', '').split('.')[-1]
    return filename


def main():
    parser = argparse.ArgumentParser(
        description="Count loops from mustache loop detection files and generate summary statistics"
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input loop TSV files (one per resolution)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output TSV file with loop counts summary'
    )
    parser.add_argument(
        '--library',
        default=None,
        help='Library name (for output)'
    )
    parser.add_argument(
        '--filter',
        default=None,
        help='Filter name (for output)'
    )
    
    args = parser.parse_args()
    
    # Count loops in each file
    loop_counts = []
    total_count = 0
    
    for loop_file in args.input:
        if not os.path.exists(loop_file):
            print(f"Warning: File not found: {loop_file}", file=sys.stderr)
            continue
        
        resolution = extract_resolution_from_filename(loop_file)
        count = count_loops_in_file(loop_file)
        loop_counts.append((resolution, count))
        total_count += count
        print(f"Resolution {resolution}: {count} loops", file=sys.stderr)
    
    # Sort by resolution (numeric if possible)
    def sort_key(x):
        try:
            return int(x[0])
        except ValueError:
            return x[0]
    
    loop_counts.sort(key=sort_key)
    
    # Write output
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(args.output, 'w') as f:
        f.write("resolution\tloop_count\n")
        for resolution, count in loop_counts:
            f.write(f"{resolution}\t{count}\n")
        f.write(f"total\t{total_count}\n")
    
    print(f"Loop counting complete. Results written to {args.output}")
    print(f"Total loops across all resolutions: {total_count}")


if __name__ == '__main__':
    main()


