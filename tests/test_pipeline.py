"""Exercise input precedence, exact count BigWigs, and end-to-end rescue."""

import csv
import json
import os
import shutil
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyBigWig
import pysam

from pooled_peaks.calibration import SHIFTS_KB, empty_counts, resolve_shifts, summarize
from pooled_peaks.cli import main
from pooled_peaks.consensus import run_consensus
from pooled_peaks.evidence import support_paths
from pooled_peaks.inputs import SAMPLE_COLUMNS, Lab, Sample, read_genome, read_samples, resolve_min_labs
from pooled_peaks.midpoints import cache_sample, count_bam, count_bigwig, export_bigwig
from pooled_peaks.pooling import direct_cores, run_pooling


MSPC = shutil.which(os.environ.get("MEPA_MSPC", "mspc"))


def make_bam(path, lengths, fragments):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": size} for name, size in lengths.items()]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as stream:
        for i, (tid, center, flag, mapq) in enumerate(fragments):
            for is_first in (True, False):
                read = pysam.AlignedSegment()
                read.query_name = f"fragment{i}"
                read.query_sequence = "A" * 50
                read.flag = flag if is_first else 147
                read.reference_id = read.next_reference_id = tid
                read.reference_start, read.next_reference_start = (
                    (center - 50, center) if is_first else (center, center - 50)
                )
                read.mapping_quality = mapq
                read.cigar = ((0, 50),)
                read.template_length = 100 if is_first else -100
                stream.write(read)


class MidpointTests(unittest.TestCase):
    def test_bam_bigwig_round_trip_and_precedence(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            lengths = {"chrTest": 2079}
            bam, bw = root / "test.bam", root / "test.bw"
            make_bam(bam, lengths, [(0, 125, 99, 0), (0, 400, 99, 20), (0, 2010, 99, 30),
                                    (0, 500, 99 | 0x400, 20), (0, 650, 99 & ~0x2, 20)])
            counts, records, depth, mapq0 = count_bam(bam, lengths)
            self.assertEqual((records, depth, mapq0), (10, 3, 1))
            self.assertEqual(np.flatnonzero(counts).tolist(), [0, 1, 8])
            self.assertEqual(export_bigwig(bam, lengths, bw), (3, 1))
            reread, total = count_bigwig(bw, lengths)
            self.assertEqual(total, depth)
            np.testing.assert_array_equal(reread, counts)
            sample = Sample("sample1", bam, root / "unselected.bw", None, None)
            report = cache_sample((sample, lengths, root))
            self.assertEqual(report["source_type"], "bam")
            self.assertEqual(report["mapq0_read1"], 1)

    def test_reject_cpm_like_bigwig(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "coverage.bw"
            with pyBigWig.open(str(path), "w") as stream:
                stream.addHeader([("chrTest", 1000)])
                stream.addEntries(["chrTest"], [0], ends=[250], values=[1.25])
            with self.assertRaisesRegex(ValueError, "raw integer midpoint count"):
                count_bigwig(path, {"chrTest": 1000})


class MultiLabTests(unittest.TestCase):
    def test_support_minimum_and_shift_counts(self):
        for count, minimum in ((2, 2), (3, 2), (4, 3), (5, 3), (6, 4), (12, 7)):
            self.assertEqual(resolve_min_labs(count, None), minimum)
            self.assertEqual(resolve_min_labs(count, count), count)
            shifts = resolve_shifts(count)
            self.assertEqual(len(set(shifts)), 6)
            for trial in shifts:
                self.assertEqual(len(set(trial)), count - 1)
                self.assertNotIn(0, trial)
        self.assertEqual(resolve_shifts(3), SHIFTS_KB)
        self.assertEqual(resolve_shifts(2), tuple((left,) for left, _ in SHIFTS_KB))
        with self.assertRaisesRegex(ValueError, "--min-labs"):
            resolve_min_labs(4, 5)
        with self.assertRaisesRegex(ValueError, "3 distinct nonzero offsets"):
            resolve_shifts(4, SHIFTS_KB)

    def test_shared_short_overlap_depends_on_required_lab_count(self):
        tracks = [
            {"chrTest": [(0, 900)]},
            {"chrTest": [(550, 1500)]},
            {"chrTest": [(550, 900)]},
            {"chrTest": [(100, 1400)]},
        ]
        self.assertEqual(direct_cores(tracks, "chrTest", 2), ([(100, 1400)], [(550, 900)]))
        self.assertEqual(direct_cores(tracks, "chrTest", 3), ([], [(550, 900)]))
        self.assertEqual(direct_cores(tracks[:2], "chrTest", 2), ([], []))

    def test_shift_adjustment_by_number_of_calling_labs(self):
        labs = [Lab(name, (), ()) for name in "ABCD"]
        counts = empty_counts(resolve_shifts(4), 4, 3)
        counts[0, 0, 0, 0, 0] = [(2, 20), (3, 10)]
        counts[0, 0, 0, 0, 1:] = [(1, 10), (1, 2)]
        first = next(summarize(counts, labs, "calibration"))
        self.assertEqual(first["observed_tiles"], 5)
        self.assertEqual(first["shifted_mean_tiles"], 7.0)
        self.assertEqual(first["shifted_observed_ratio"], 1.4)

    def test_weak_support_needs_one_extra_lab(self):
        strong, extra_weak = support_paths(np.array([1, 1, 2, 2, 3]), np.array([1, 1, 2, 2, 2]),
                                           np.array([1.0, 2.0, 1.0, 2.0, 1.0]), 2.0)
        self.assertEqual(strong.tolist(), [False, True, False, True, False])
        self.assertEqual(extra_weak.tolist(), [False, False, False, False, True])

    def test_two_four_and_six_lab_cli_rescue_and_override(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            genome = root / "genome.tsv"
            genome.write_text("chrom\tlength\tsplit\nchrFit\t5000000\tcalibration\n")
            sheet = root / "samples.tsv"
            rows = []
            shared_by_three, shared_by_all = 1_700_000, 1_800_000
            single_call, double_call = 2_000_000, 2_200_000
            for lab in "ABCDEF":
                called = [shared_by_all]
                if lab in "ABC":
                    called.append(shared_by_three)
                if lab == "A":
                    called.append(single_call)
                if lab in "AB":
                    called.append(double_call)
                near = [single_call, double_call] if lab in "ABCE" else []
                primary, permissive = root / f"{lab}.broadPeak", root / f"{lab}.permissive.bed"
                primary.write_text("".join(
                    f"chrFit\t{start}\t{start + 3000}\t{lab}_{start}\t900\t.\t9\t9\t9\n"
                    for start in sorted(called)
                ))
                permissive.write_text("".join(
                    f"chrFit\t{start}\t{start + 3000}\t{lab}_{start}\n" for start in sorted(set(called + near))
                ))
                midpoint = root / f"{lab}.bw"
                counts = {1_500_000: 1}
                for start in near:
                    profile = (1, 2, 3, 7, 9, 13, 19, 14, 10, 5, 3, 2, 1)
                    counts.update({start + (i - 2) * 250: value for i, value in enumerate(profile)})
                starts = sorted(counts)
                with pyBigWig.open(str(midpoint), "w") as stream:
                    stream.addHeader([("chrFit", 5_000_000)])
                    stream.addEntries(["chrFit"] * len(starts), starts,
                                      ends=[start + 250 for start in starts],
                                      values=[float(counts[start]) for start in starts])
                rows.append([lab, lab, "target", primary, permissive, "", midpoint])

            def run(count):
                with sheet.open("w", newline="") as stream:
                    writer = csv.writer(stream, delimiter="\t")
                    writer.writerow(SAMPLE_COLUMNS)
                    writer.writerows(rows[:count])
                output = root / f"labs{count}"
                with redirect_stdout(StringIO()):
                    code = main(["run", "--samples", str(sheet), "--genome", str(genome),
                                 "--workers", "1", "--outdir", str(output)])
                self.assertEqual(code, 0)
                return output

            two = run(2)
            six = run(6)
            four = run(4)
            for run_root, minimum in ((two, 2), (four, 3), (six, 4)):
                output = run_root / "step4_pool"
                thresholds = json.loads((output / "thresholds.json").read_text())
                self.assertEqual(thresholds["status"], "calibrated")
                self.assertEqual(thresholds["min_labs"], minimum)
                self.assertEqual(thresholds["calibration"]["shifted_observed_ratio"], 0)
                with (output / "pooled_peak_evidence.tsv").open() as stream:
                    regions = list(csv.DictReader(stream, delimiter="\t"))
                for locus in (single_call, double_call) if minimum > 2 else (single_call,):
                    rescued = next(row for row in regions if int(row["start"]) <= locus < int(row["end"]))
                    self.assertGreater(int(rescued["added_signal_seed_bp"]), 0)

            strict = root / "unanimous"
            with redirect_stdout(StringIO()):
                code = main(["step4", "--samples", str(sheet), "--genome", str(genome),
                             "--step3-dir", str(four / "step3_mspc"), "--outdir", str(strict),
                             "--workers", "1", "--min-labs", "4"])
            self.assertEqual(code, 0)
            summary = json.loads((strict / "pooling_summary.json").read_text())
            self.assertEqual((summary["lab_count"], summary["min_labs"], summary["peaks"]), (4, 4, 1))
            self.assertTrue((strict / "pooled_peaks.bed").read_text().startswith(f"chrFit\t{shared_by_all}\t"))


@unittest.skipUnless(MSPC, "MSPC v6.1.0 is required for this integration test")
class PipelineTests(unittest.TestCase):
    def test_mspc_auto_cutoff_and_same_output_from_bigwig(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome.tsv").write_text("chrom\tlength\tsplit\n"
                                             "chrFit\t3000000\tcalibration\n"
                                             "chrCheck\t3000000\tholdout\n")
            lengths, splits = read_genome(root / "genome.tsv")
            layout = [("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2"), ("C", "C1"), ("A", "IgG")]
            bam_rows, bw_rows = [], []
            for lab, name in layout:
                control = name == "IgG"
                bam, bw = root / f"{name}.bam", root / f"{name}.bw"
                primary = root / f"{name}.broadPeak"
                permissive = root / f"{name}.permissive.broadPeak"
                peak_start = {"A": 1_300_000, "B": 1_301_000, "C": 1_400_000}[lab]
                fragment_rows, peak_rows, support_rows = [], [], []
                for tid, chrom in enumerate(lengths):
                    fragment_rows.append((tid, 1_000_000, 99, 0))
                    if control:
                        fragment_rows.extend((tid, 1_010_000 + 250 * i, 99, 0) for i in range(40))
                        continue
                    peak_rows.append(f"{chrom}\t{peak_start}\t{peak_start + 3000}\t{name}_{tid}\t"
                                     "900\t.\t9\t9\t9\n")
                    if lab == "A":
                        support_rows.append(f"{chrom}\t1399800\t1403500\t{name}_nearby_{tid}\n")
                    if lab in ("A", "C"):
                        heights = (1, 2, 3, 5, 9, 14, 17, 13, 8, 5, 3, 2, 1)
                        for index, number in enumerate(heights):
                            fragment_rows.extend((tid, 1_399_875 + 250 * index, 99, 0) for _ in range(number))
                make_bam(bam, lengths, fragment_rows)
                export_bigwig(bam, lengths, bw)
                role = "igg" if control else "target"
                if not control:
                    primary.write_text("".join(peak_rows))
                    permissive.write_text("".join(support_rows))
                metadata = [lab, name, role, "" if control else primary, permissive if support_rows else ""]
                bam_rows.append([*metadata, bam, bw])
                bw_rows.append([*metadata, "", bw])

            for sheet, rows in (("bams.tsv", bam_rows), ("bigwigs.tsv", bw_rows)):
                with (root / sheet).open("w", newline="") as stream:
                    writer = csv.writer(stream, delimiter="\t")
                    writer.writerow(SAMPLE_COLUMNS)
                    writer.writerows(rows)

            bam_labs = read_samples(root / "bams.tsv")
            step3 = root / "step3"
            run_consensus(bam_labs, lengths, step3, str(Path(MSPC).resolve()), 1)
            self.assertIn("chrFit\t1300000\t1303000", (step3 / "peak_beds/A.MSPC_consensus.bed").read_text())
            results = []
            for sheet in ("bams.tsv", "bigwigs.tsv"):
                out = root / f"step4_{sheet}"
                report = run_pooling(read_samples(root / sheet), lengths, splits,
                                     step3, out, 1, 0.10, SHIFTS_KB)
                self.assertEqual(report["status"], "calibrated")
                self.assertGreater(report["added_signal_seed_bp"], 0)
                self.assertEqual(report["direct_500_bp"], 4000)
                thresholds = json.loads((out / "thresholds.json").read_text())
                self.assertEqual(thresholds["calibration"]["shifted_observed_ratio"], 0)
                self.assertEqual(thresholds["holdout"]["shifted_observed_ratio"], 0)
                results.append(tuple((out / name).read_bytes() for name in (
                    "pooled_peaks.bed", "pooled_peak_evidence.tsv", "thresholds.json",
                    "signal_threshold_calibration.tsv",
                )))
            self.assertEqual(*results)
