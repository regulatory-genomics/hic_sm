import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "workflow"
    / "scripts"
    / "hictailor_json_to_multiqc_summary.py"
)
DEDUP_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "workflow"
    / "scripts"
    / "hictailor_json_to_dedup_stats.py"
)


def load_module(path=SCRIPT, name="hictailor_multiqc"):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class HicTailorMultiqcSummaryTests(unittest.TestCase):
    def test_required_depth_returns_pcr_reads_and_blank_total_depth(self):
        module = load_module()

        n_pcr, n_total = module.calculate_required_total_reads(
            target_u=2_000_000_000,
            complexity=3_000_000_000,
            optical_rate=0.25,
            include_total_depth=False,
        )

        expected = -3_000_000_000 * math.log(1 - 2_000_000_000 / 3_000_000_000)
        self.assertAlmostEqual(n_pcr, expected, places=3)
        self.assertEqual(n_total, "")

    def test_required_depth_returns_na_when_target_exceeds_complexity(self):
        module = load_module()

        n_pcr, n_total = module.calculate_required_total_reads(
            target_u=2_000_000_000,
            complexity=1_000_000_000,
            optical_rate=0.1,
            include_total_depth=False,
        )

        self.assertEqual(n_pcr, "NA")
        self.assertEqual(n_total, "")

    def test_estimate_unique_2b_uses_linear_log_formula_for_h8_12(self):
        module = load_module()

        estimate = module.estimate_unique_2b_without_optical(
            initial_raw_pair=43_617_627,
            initial_vp_unique=12_330_361,
        )

        expected = math.expm1(
            11.014228
            - 2.730083540102 * math.log1p(43_617_627)
            + 3.309305076052 * math.log1p(12_330_361)
        )
        self.assertEqual(module.format_count(estimate), "24700465")
        self.assertEqual(module.format_count(estimate), module.format_count(expected))

    def test_effective_discard_rate_counts_stitched_uncut_when_available(self):
        module = load_module()

        stats = {
            "total_pairs_read": 1000,
            "stitched_uncut": 400,
            "unstitched_discarded": 25,
            "rates": {"total_discard_rate": 0.025},
        }

        self.assertAlmostEqual(module.effective_discard_rate(stats), 0.425)

    def test_build_summary_row_uses_hictailor_json_pairs_dist_and_loop_inputs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            stats_json = tmpdir / "sample.json"
            pairs = tmpdir / "sample.pairs"
            dist = tmpdir / "dist.csv"
            loops = tmpdir / "loops.tsv"

            stats_json.write_text(
                json.dumps(
                    {
                        "total_pairs_read": 1000,
                        "valid_pairs": 400,
                        "trans": 100,
                        "hic_dangling_end": 20,
                        "hic_self_circle_frag": 10,
                        "stitched": 500,
                        "rates": {
                            "total_discard_rate": 0.123,
                            "align_percentage": 0.8,
                        },
                        "library_complexity": {
                            "total_mapped_for_complexity": 500,
                            "total_dups_for_complexity": 50,
                            "complexity_dups_by_tile_median": 3_000_000_000,
                            "estimated_optical_duplicates": 25,
                        },
                    }
                )
            )
            pairs.write_text(
                "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n"
                "r1 chr1 1000 chr1 25000 + - UU 60 60\n"
                "r2 chr1 1000 chr1 15000 + - UU 60 60\n"
                "r3 chr1 1000 chr2 2000 + - UU 60 60\n"
            )
            dist.write_text("region,slope,mse,n_points\nchr1p,-1.0,0.3,10\nALL_REGIONS,-1.1,0.42,99\n")
            loops.write_text("resolution\tloop_count\n5000\t7\ntotal\t11\n")

            row = module.build_summary_row(
                library="sample",
                stats_json=stats_json,
                pairs_path=pairs,
                dist_path=dist,
                loop_path=loops,
                target_unique=2_000_000_000,
                include_total_depth=False,
            )

        self.assertEqual(row["sample"], "sample")
        self.assertEqual(row["raw_pair"], 1000)
        self.assertEqual(row["VP_unique"], 400)
        self.assertAlmostEqual(row["VP_rate"], 0.4)
        self.assertAlmostEqual(row["Discard rate"], 0.123)
        self.assertEqual(row["VP_cis>20K"], 1)
        self.assertAlmostEqual(row["VP_cis_20k/valid_pair"], 1 / 400)
        self.assertAlmostEqual(row["trans/VP_unique"], 0.25)
        self.assertAlmostEqual(row["cis/trans ratio"], 3.0)
        self.assertAlmostEqual(row["Dangling%"], 0.02)
        self.assertAlmostEqual(row["Self_ligation%"], 0.01)
        self.assertAlmostEqual(row["linker_rate"], 0.5)
        self.assertAlmostEqual(row["Map%"], 0.8)
        self.assertAlmostEqual(row["Dup%"], 0.1)
        self.assertAlmostEqual(row["dist_mse_log"], 0.42)
        self.assertEqual(row["Loop_count"], 11)
        self.assertEqual(row["estimate_sequence_depth"], "")
        expected_unique_2b = module.estimate_unique_2b_without_optical(1000, 400)
        self.assertEqual(
            row["estimate_unique_2B_without_optical"],
            module.format_count(expected_unique_2b),
        )

    def test_dedup_stats_bridge_uses_hictailor_library_complexity(self):
        module = load_module(DEDUP_SCRIPT, "hictailor_dedup_stats")
        data = {
            "total_pairs_read": 1000,
            "valid_pairs": 400,
            "pcr_dup_count": 12,
            "optical_dup_count": 3,
            "exact_dedup_hits": 5,
            "library_complexity": {
                "total_mapped_for_complexity": 500,
                "total_dups_for_complexity": 50,
                "complexity_naive": 12345.5,
                "dups_by_tile_median": 44,
                "estimated_optical_duplicates": 6,
                "complexity_dups_by_tile_median": 23456.5,
                "tile_pairs_observed": 7,
                "tile_parse_failed": 8,
            },
        }

        rows = dict(module.build_rows(data))

        self.assertEqual(rows["total"], 1000)
        self.assertEqual(rows["summary/total_mapped_for_complexity"], 500)
        self.assertEqual(rows["summary/total_dups_for_complexity"], 50)
        self.assertAlmostEqual(rows["summary/frac_dups"], 0.1)
        self.assertAlmostEqual(rows["summary/optical_rate"], 6 / 500)
        self.assertEqual(rows["summary/complexity_naive"], 12345.5)
        self.assertEqual(rows["summary/dups_by_tile_median"], 44)
        self.assertEqual(rows["summary/complexity_dups_by_tile_median"], 23456.5)
        self.assertEqual(rows["summary/tile_pairs_observed"], 7)
        self.assertEqual(rows["summary/tile_parse_failed"], 8)


if __name__ == "__main__":
    unittest.main()
