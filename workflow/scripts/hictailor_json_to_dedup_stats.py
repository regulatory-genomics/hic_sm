#!/usr/bin/env python3
"""
Convert hic-tailor JSON stats into a minimal pairtools-like dedup.stats TSV.

This is a compatibility bridge so existing hic_sm downstream rules
(calculate_complexity, MultiQC inputs, report links) keep working while
hic-tailor owns global sorting + deduplication.
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="hic-tailor JSON stats")
    parser.add_argument("-o", "--output", required=True, help="output dedup.stats TSV")
    args = parser.parse_args()

    with open(args.input) as f:
        data = json.load(f)

    total_pairs_read = int(data.get("total_pairs_read", 0))
    valid_pairs = int(data.get("valid_pairs", data.get("hic_valid", 0)))
    pcr_dup_count = int(data.get("pcr_dup_count", 0))
    optical_dup_count = int(data.get("optical_dup_count", 0))
    exact_dedup_hits = int(data.get("exact_dedup_hits", 0))

    # Compatibility approximations:
    # - total: keep pairtools-like total as total read pairs entering pipeline
    # - optical rate: derive from valid pairs when available
    # - complexity_dups_by_tile_median: fallback proxy so calculate_complexity.py can run
    #   We use deduplicated valid pairs as a conservative lower-bound proxy.
    optical_rate = (optical_dup_count / valid_pairs) if valid_pairs > 0 else 0.0
    complexity_proxy = max(valid_pairs, 1)

    rows = [
        ("total", total_pairs_read),
        ("summary/frac_dups", (pcr_dup_count / valid_pairs) if valid_pairs > 0 else 0.0),
        ("summary/optical_rate", optical_rate),
        ("summary/complexity_dups_by_tile_median", complexity_proxy),
        ("summary/hictailor_valid_pairs", valid_pairs),
        ("summary/hictailor_pcr_dup_count", pcr_dup_count),
        ("summary/hictailor_optical_dup_count", optical_dup_count),
        ("summary/hictailor_exact_dedup_hits", exact_dedup_hits),
    ]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        for key, value in rows:
            f.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    main()
