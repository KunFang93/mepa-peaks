"""Exercise input precedence, exact count BigWigs, and end-to-end rescue."""

import csv
import json
import os
import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyBigWig
import pysam

from pooled_peaks.calibration import SHIFTS_KB
from pooled_peaks.consensus import run_consensus
from pooled_peaks.inputs import SAMPLE_COLUMNS, Sample, read_genome, read_samples
from pooled_peaks.midpoints import cache_sample, count_bam, count_bigwig, export_bigwig
from pooled_peaks.pooling import run_pooling


MSPC = shutil.which(os.environ.get("CONFIDENT_PEAKS_MSPC", "mspc"))


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
