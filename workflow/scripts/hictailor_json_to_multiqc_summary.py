#!/usr/bin/env python3
from __future__ import annotations
"""Build a compact hic-tailor summary table for MultiQC custom data.

The script consumes hic-tailor's JSON report plus small downstream summaries and
emits one TSV row per library. It keeps workflow logic out of Snakemake rules and
centralizes the metric definitions used by the final MultiQC table.
"""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any


FIELDS = [
    "sample",
    "raw_pair",
    "VP_rate",
    "VP_unique",
    "Discard rate",
    "VP_cis>20K",
    "Map%",
    "Dup%",
    "VP_cis_20k/valid_pair",
    "trans/VP_unique",
    "cis/trans ratio",
    "Dangling%",
    "Self_ligation%",
    "linker_rate",
    "estimate_unique_2B_without_optical",
    "estimate_sequence_depth",
    "dist_mse_log",
    "Loop_count",
]


def safe_div(numerator: float, denominator: float) -> float:
    return 0.0 if denominator == 0 else numerator / denominator


def format_count(value: float | str) -> str:
    if value == "NA":
        return "NA"
    return str(int(round(float(value))))


def effective_discard_rate(stats: dict[str, Any]) -> float | str:
    """Return cut-mode discard rate when stitched-uncut counts are present.

    hic-tailor CLI JSON keeps the legacy Micro-C-friendly value in
    rates.total_discard_rate. This workflow runs restriction-enzyme cut-mode
    Hi-C, where stitched uncut reads are not emitted as valid pairs, so the
    report should include them in Discard rate whenever raw counters exist.
    """
    total_pairs = int(stats.get("total_pairs_read", 0) or 0)
    if total_pairs and "stitched_uncut" in stats and "unstitched_discarded" in stats:
        return (
            int(stats.get("stitched_uncut", 0) or 0)
            + int(stats.get("unstitched_discarded", 0) or 0)
        ) / total_pairs
    return (stats.get("rates", {}) or {}).get("total_discard_rate", "")


def load_json(path: Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        return json.load(handle)


def open_text_auto(path: Path):
    path = Path(path)
    with path.open("rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return path.open()


def estimate_unique_2b_without_optical(
    initial_raw_pair: float,
    initial_vp_unique: float,
) -> float | str:
    """Predict unique reads at 2B sequences from initial raw/valid pairs.

    This is a forward estimate requested for the MultiQC column
    `estimate_unique_2B_without_optical`; it does not use the old inverse
    library-complexity depth calculation.
    """
    raw_pair = float(initial_raw_pair or 0)
    vp_unique = float(initial_vp_unique or 0)
    if raw_pair <= 0 or vp_unique <= 0:
        return "NA"
    estimate = math.expm1(
        11.014228
        - 2.730083540102 * math.log1p(raw_pair)
        + 3.309305076052 * math.log1p(vp_unique)
    )
    if not math.isfinite(estimate) or estimate < 0:
        return "NA"
    return estimate


def calculate_required_total_reads(
    target_u: float,
    complexity: float,
    optical_rate: float,
    include_total_depth: bool = False,
) -> tuple[float | str, float | str]:
    """Return required effective PCR reads and optional total reads.

    Formula: U = C * (1 - exp(-N/C)); therefore N = -C * ln(1 - U/C).
    If target_u >= C, the target unique count is mathematically unreachable and
    the effective read estimate is reported as NA. The current workflow leaves
    estimate_sequence_depth blank by request, so include_total_depth defaults to
    False.
    """
    if complexity <= 0:
        return "NA", "" if not include_total_depth else "NA"
    if target_u >= complexity:
        return "NA", "" if not include_total_depth else "NA"
    if optical_rate < 0 or optical_rate >= 1:
        return "NA", "" if not include_total_depth else "NA"

    n_pcr = -complexity * math.log(1 - target_u / complexity)
    if not include_total_depth:
        return n_pcr, ""
    return n_pcr, n_pcr / (1 - optical_rate)


def count_cis_pairs_over_distance(pairs_path: Path, min_distance: int = 20_000) -> int:
    """Count cis pairs with genomic distance greater than min_distance.

    The workflow writes 4DN-style pairs where columns 2-5 are chrom1, pos1,
    chrom2, pos2. Extra columns such as parent_readID are ignored.
    """
    count = 0
    with open_text_auto(Path(pairs_path)) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 5:
                continue
            chrom1, chrom2 = fields[1], fields[3]
            if chrom1 != chrom2:
                continue
            try:
                pos1, pos2 = int(fields[2]), int(fields[4])
            except ValueError:
                continue
            if abs(pos2 - pos1) > min_distance:
                count += 1
    return count


def read_all_regions_mse(path: Path) -> float | str:
    if not Path(path).exists():
        return ""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("region") == "ALL_REGIONS":
                value = row.get("mse", "")
                return "" if value in {"", "nan", "NaN"} else float(value)
    return ""


def read_total_loop_count(path: Path) -> int | str:
    if not Path(path).exists():
        return ""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("resolution") == "total":
                return int(float(row.get("loop_count", 0)))
    return ""


def build_summary_row(
    library: str,
    stats_json: Path,
    pairs_path: Path,
    dist_path: Path,
    loop_path: Path,
    target_unique: float = 2_000_000_000,
    include_total_depth: bool = False,
) -> dict[str, Any]:
    data = load_json(Path(stats_json))
    rates = data.get("rates", {}) or {}
    complexity = data.get("library_complexity", {}) or {}

    raw_pair = int(data.get("total_pairs_read", 0) or 0)
    valid_pairs = int(data.get("valid_pairs", data.get("hic_valid", 0)) or 0)
    trans = int(data.get("trans", 0) or 0)
    cis_pairs = max(valid_pairs - trans, 0)
    cis_over_20k = count_cis_pairs_over_distance(Path(pairs_path), 20_000)

    total_mapped = float(complexity.get("total_mapped_for_complexity", 0) or 0)
    total_dups = float(complexity.get("total_dups_for_complexity", 0) or 0)
    estimate_unique_2b = estimate_unique_2b_without_optical(raw_pair, valid_pairs)
    estimate_total = ""

    return {
        "sample": library,
        "raw_pair": raw_pair,
        "VP_rate": safe_div(valid_pairs, raw_pair),
        "VP_unique": valid_pairs,
        "Discard rate": effective_discard_rate(data),
        "VP_cis>20K": cis_over_20k,
        "Map%": rates.get("align_percentage", rates.get("mapping_rate", "")),
        "Dup%": safe_div(total_dups, total_mapped),
        "VP_cis_20k/valid_pair": safe_div(cis_over_20k, valid_pairs),
        "trans/VP_unique": safe_div(trans, valid_pairs),
        "cis/trans ratio": "" if trans == 0 else cis_pairs / trans,
        "Dangling%": safe_div(int(data.get("hic_dangling_end", 0) or 0), raw_pair),
        "Self_ligation%": safe_div(int(data.get("hic_self_circle_frag", 0) or 0), raw_pair),
        "linker_rate": safe_div(int(data.get("stitched", 0) or 0), raw_pair),
        "estimate_unique_2B_without_optical": format_count(estimate_unique_2b),
        "estimate_sequence_depth": estimate_total,
        "dist_mse_log": read_all_regions_mse(Path(dist_path)),
        "Loop_count": read_total_loop_count(Path(loop_path)),
    }


def write_summary(row: dict[str, Any], output: Path) -> None:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create hic-tailor MultiQC custom summary TSV")
    parser.add_argument("--library", required=True)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--dist", required=True, type=Path)
    parser.add_argument("--loops", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-unique", type=float, default=2_000_000_000)
    parser.add_argument(
        "--include-total-depth",
        action="store_true",
        help="Also fill estimate_sequence_depth with optical-rate-corrected total reads",
    )
    args = parser.parse_args()

    row = build_summary_row(
        library=args.library,
        stats_json=args.json,
        pairs_path=args.pairs,
        dist_path=args.dist,
        loop_path=args.loops,
        target_unique=args.target_unique,
        include_total_depth=args.include_total_depth,
    )
    write_summary(row, args.output)


if __name__ == "__main__":
    main()
