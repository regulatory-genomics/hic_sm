#!/usr/bin/env python3
"""
Calculate library complexity (expected unique observed reads) from pairtools dedup statistics.

Formula:
  N_optical = N_optical_rate * N_total
  N_pcr = N_total - N_optical
  U_observed = C * (1 - exp(-N_pcr/C))

Where:
  - N_total: total number of reads
  - C: library complexity (from summary/complexity_dups_by_tile_median)
  - N_optical_rate: optical duplicate rate (fraction)
  - U_observed: expected number of observed unique reads
"""

import argparse
import sys
import math
import os
from pathlib import Path


def parse_stats_file(stats_file):
    """Parse pairtools dedup stats file (tab-separated key-value pairs)."""
    stats = {}
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '\t' not in line:
                    continue
                key, value = line.split('\t', 1)
                try:
                    # Try to convert to float, fall back to string
                    stats[key] = float(value)
                except ValueError:
                    stats[key] = value
    except Exception as e:
        print(f"Error reading stats file {stats_file}: {e}", file=sys.stderr)
        sys.exit(1)
    return stats


def calculate_unique_reads_optical(N_total, C, N_optical_rate):
    """
    Calculate expected number of observed unique reads.
    
    Args:
        N_total: Total number of reads
        C: Library complexity
        N_optical_rate: Optical duplicate rate (fraction, 0-1)
    
    Returns:
        U_observed: Expected number of observed unique reads
    """
    # Input validation
    if not all(isinstance(x, (int, float)) for x in [N_total, C, N_optical_rate]):
        raise ValueError("Inputs N_total, C, and N_optical_rate must be numeric.")
    
    if C <= 0:
        raise ValueError("C (library complexity) must be a positive number.")
    
    if N_total < 0:
        print(f"Warning: N_total ({N_total}) is negative. The result may be nonsensical.", file=sys.stderr)
    
    # Calculate N_optical
    N_optical = N_optical_rate * N_total
    
    if N_optical > N_total:
        raise ValueError(f"N_optical ({N_optical}) cannot be greater than N_total ({N_total}).")
    
    # Calculate N_pcr (the "effective" N for the complexity formula)
    # This is the total number of reads *minus* the technical artifacts.
    N_pcr = N_total - N_optical
    
    # Calculate the ratio -N_pcr/C
    ratio = -N_pcr / C
    
    # Implement the formula: U_observed = C * (1 - exp(-N_pcr/C))
    # This U is the expected number of *observed unique reads*.
    U_observed = C * (1 - math.exp(ratio))
    
    return U_observed, N_optical, N_pcr


def main():
    parser = argparse.ArgumentParser(
        description="Calculate library complexity from pairtools dedup statistics"
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input dedup.stats file from pairtools'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output TSV file with complexity metrics'
    )
    parser.add_argument(
        '--optical-rate',
        type=float,
        default=None,
        help='Optical duplicate rate (fraction, 0-1). If not provided, will try to estimate from stats or use 0'
    )
    parser.add_argument(
        '--library',
        default=None,
        help='Library name (will be extracted from input filename if not provided)'
    )
    
    args = parser.parse_args()
    
    # Parse stats file
    stats = parse_stats_file(args.input)
    
    # Extract required values
    if 'total' not in stats:
        print("Error: 'total' field not found in stats file", file=sys.stderr)
        sys.exit(1)
    N_total = stats['total']
    
    if 'summary/complexity_dups_by_tile_median' not in stats:
        print("Error: 'summary/complexity_dups_by_tile_median' field not found in stats file", file=sys.stderr)
        sys.exit(1)
    C = stats['summary/complexity_dups_by_tile_median']
    
    # Determine optical rate
    if args.optical_rate is not None:
        N_optical_rate = args.optical_rate
    elif 'summary/optical_rate' in stats:
        N_optical_rate = stats['summary/optical_rate']
    elif 'dups_by_tile_median' in stats and N_total > 0:
        # Estimate from tile duplicates if available
        N_optical_rate = stats['dups_by_tile_median'] / N_total
        print(f"Warning: Using estimated optical rate from tile duplicates: {N_optical_rate:.6f}", file=sys.stderr)
    else:
        # Default to 0 if not available
        N_optical_rate = 0.0
        print(f"Warning: Optical rate not found, using default: {N_optical_rate}", file=sys.stderr)
    
    # Calculate complexity
    try:
        U_observed, N_optical, N_pcr = calculate_unique_reads_optical(N_total, C, N_optical_rate)
    except Exception as e:
        print(f"Error calculating complexity: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Get library name
    if args.library is None:
        library = Path(args.input).stem.replace('.dedup.stats', '').replace('.dedup', '')
    else:
        library = args.library
    
    # Write output
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(args.output, 'w') as f:
        f.write("library\tN_total\tC\tN_optical_rate\tN_optical\tN_pcr\tU_observed\n")
        f.write(f"{library}\t{N_total:.0f}\t{C:.2f}\t{N_optical_rate:.6f}\t{N_optical:.0f}\t{N_pcr:.0f}\t{U_observed:.2f}\n")
    
    print(f"Complexity calculation complete. Results written to {args.output}")
    print(f"Library: {library}")
    print(f"N_total: {N_total:.0f}")
    print(f"C (complexity): {C:.2f}")
    print(f"N_optical_rate: {N_optical_rate:.6f}")
    print(f"N_optical: {N_optical:.0f}")
    print(f"N_pcr: {N_pcr:.0f}")
    print(f"U_observed: {U_observed:.2f}")


if __name__ == '__main__':
    main()


