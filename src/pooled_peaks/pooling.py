"""Build direct cores, rescue single-lab calls, and grow supported boundaries."""

from collections import defaultdict
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from .calibration import METRIC_COLUMNS, calibrate_chrom, empty_counts, select_thresholds, summarize
from .consensus import checked_consensus
from .evidence import (
    MAX_LAG_BINS, PERMISSIVE_PAD, SCORE_HALF, aligned_shape, build_features, flags,
    normalize_profiles, sustained,
)
from .inputs import BIN, read_peaks, sha256, stage_directory, write_json, write_tsv
from .midpoints import MidpointCounts


MIN_DIRECT = 500
MIN_SHORT_OVERLAP = 250
SHORTER_FRACTION = 0.5
EVIDENCE_FIELDS = ("chrom", "start", "end", "peak", "bp", "direct_500_bp", "added_short_direct_bp",
                   "added_signal_seed_bp", "added_shoulder_bp", "strong_peer_seed_bp", "both_weak_seed_bp")


def union(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def overlaps(first, second):
    i = j = 0
    while i < len(first) and j < len(second):
        a, b = first[i]
        c, d = second[j]
        start, end = max(a, c), min(b, d)
        if start < end:
            yield start, end, b - a, d - c
        i += b <= d
        j += d <= b


def direct_cores(tracks, chrom):
    all_overlaps, short = [], []
    for first, second in combinations(tracks, 2):
        for a, b, w1, w2 in overlaps(first.get(chrom, []), second.get(chrom, [])):
            width = b - a
            all_overlaps.append((a, b))
            if MIN_SHORT_OVERLAP <= width < MIN_DIRECT and width >= SHORTER_FRACTION * min(w1, w2):
                short.append((a, b))
    long = [(a, b) for a, b in union(all_overlaps) if b - a >= MIN_DIRECT]
    return long, union(short)


def spans_by_chrom(peaks, disjoint=False):
    tracks = defaultdict(list)
    for peak in peaks:
        intervals = tracks[peak.chrom]
        if disjoint and intervals and intervals[-1][1] > peak.start:
            raise ValueError(f"Lab-level peaks overlap at {peak.chrom}:{peak.start}-{peak.end}")
        intervals.append((peak.start, peak.end))
    return dict(tracks)


def tile_runs(mask, length):
    corners = np.flatnonzero(np.diff(np.r_[False, mask, False]))
    return [(int(a) * BIN, min(int(b) * BIN, length)) for a, b in zip(corners[::2], corners[1::2])]


def signal_seeds(features, original, thresholds):
    scores, peer_scores, profiles, called = features
    n = scores.shape[1]
    weak, strong, min_corr = (thresholds[key] for key in ("weak", "strong", "min_corr"))
    eligible = (called.sum(axis=0) == 1) & ~flags(original, n, center_only=False)
    strong_hits, weak_hits = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    for i, anchor_calls in enumerate(called):
        indices = np.flatnonzero(eligible & anchor_calls & (scores[i] >= weak))
        if not len(indices):
            continue
        anchor = normalize_profiles(profiles[i], indices)
        peers = []
        for j in range(len(called)):
            if j != i:
                match = aligned_shape(anchor, profiles[j], indices) >= min_corr
                peers.append(np.where(match, peer_scores[j, indices], -np.inf))
        strong_peer = np.maximum(*peers) >= strong
        both_weak = np.minimum(*peers) >= weak
        passing = np.zeros(n, dtype=bool)
        passing[indices] = strong_peer | both_weak
        passing = sustained(passing)[indices]
        strong_hits[indices] |= passing & strong_peer
        weak_hits[indices] |= passing & both_weak
    return strong_hits, weak_hits


class Footprint:
    def __init__(self, spans):
        arr = np.asarray(spans, dtype=np.int64).reshape(-1, 2)
        self.starts = np.r_[arr[:, 0], np.iinfo(np.int64).max]
        self.ends = arr[:, 1]
        self.prefix = np.r_[0, np.cumsum(arr[:, 1] - arr[:, 0])]

    def mass(self, x):
        i = np.searchsorted(self.ends, x, side="right")
        return self.prefix[i] + np.maximum(0, x - self.starts[i])

    def coverage(self, start, end):
        return int(self.mass(end) - self.mass(start))


def pool_chrom(original, short, envelope, length, features, thresholds):
    direct = union(original + short)
    seeds = final = direct
    strong_spans, weak_spans = [], []
    if thresholds is not None:
        scores, peers, _, called = features
        strong, weak = signal_seeds(features, original, thresholds)
        signal_spans = [(a, b) for a, b, _, _ in overlaps(tile_runs(strong | weak, length), envelope)
                        if b - a >= 2 * BIN]
        seeds = union(direct + signal_spans)
        passing = (peers >= thresholds["weak"]) | (called & (scores >= thresholds["weak"]))
        shoulders = sustained((passing.sum(axis=0) >= 2) & flags(envelope, scores.shape[1]))
        seed_footprint = Footprint(seeds)
        growth = [(a, b) for a, b, _, _ in overlaps(tile_runs(shoulders, length), envelope)
                  if seed_footprint.coverage(a, b) > 0]
        final = union(seeds + growth)
        strong_spans = [(a, b) for a, b, _, _ in overlaps(tile_runs(strong, length), signal_spans)]
        weak_spans = [(a, b) for a, b, _, _ in overlaps(tile_runs(weak, length), signal_spans)]
    footprints = [Footprint(spans) for spans in (original, direct, seeds, strong_spans, weak_spans)]
    for a, b in final:
        old_bp, direct_bp, seed_bp, strong_bp, weak_bp = [fp.coverage(a, b) for fp in footprints]
        yield dict(start=a, end=b, bp=b - a, direct_500_bp=old_bp,
                   added_short_direct_bp=direct_bp - old_bp, added_signal_seed_bp=seed_bp - direct_bp,
                   added_shoulder_bp=b - a - seed_bp, strong_peer_seed_bp=strong_bp, both_weak_seed_bp=weak_bp)


def read_sources(labs, lengths, step3):
    target_samples = [sample for lab in labs for sample in lab.targets]
    primary_hashes = {sample.name: sha256(sample.broadpeak) for sample in target_samples}
    hard = [spans_by_chrom(peaks, disjoint=True)
            for peaks in checked_consensus(labs, lengths, step3, primary_hashes)]
    combined = defaultdict(list)
    permissive, pins = {}, []
    for sample in target_samples:
        for chrom, spans in spans_by_chrom(read_peaks(sample.broadpeak, lengths)).items():
            combined[chrom].extend(spans)
        supporting = sample.permissive or sample.broadpeak
        permissive[sample.name] = spans_by_chrom(read_peaks(supporting, lengths))
        pins.append(dict(sample=sample.name, primary=str(sample.broadpeak), primary_sha256=primary_hashes[sample.name],
                         permissive=str(supporting), permissive_sha256=sha256(supporting)))
    return hard, {c: union(spans) for c, spans in combined.items()}, permissive, pins


def run_pooling(labs, lengths, splits, step3, destination, workers, target, shifts):
    if "calibration" not in splits.values():
        raise ValueError("Step 4 needs at least one chromosome labeled calibration in the genome TSV")
    hard, envelope, permissive, pins = read_sources(labs, lengths, step3)
    samples = [sample for lab in labs for sample in (*lab.targets, *lab.controls)]
    geometry = {chrom: direct_cores(hard, chrom) for chrom in lengths}
    with stage_directory(destination) as root:
        write_tsv(root / "peak_inputs.tsv", tuple(pins[0]), pins)
        with TemporaryDirectory(prefix=".midpoints-", dir=root) as cache_dir:
            cache = MidpointCounts(samples, lengths, Path(cache_dir), root / "count_summary.tsv", workers)
            grid_counts = {split: empty_counts(shifts) for split in ("calibration", "holdout")
                           if split in splits.values()}
            for chrom, length in lengths.items():
                split = splits[chrom]
                if split in grid_counts:
                    features = build_features(labs, cache, hard, permissive, chrom, length)
                    grid_counts[split] += calibrate_chrom(features, geometry[chrom][0], shifts)
            metrics = [row for split, counts in grid_counts.items() for row in summarize(counts, labs, split)]
            write_tsv(root / "signal_threshold_calibration.tsv", METRIC_COLUMNS, metrics)
            thresholds = select_thresholds(metrics, target)
            settings = dict(status="calibrated" if thresholds else "direct_only",
                            target_shifted_ratio=target, splits=splits, shifts_kb=shifts, bin_bp=BIN,
                            score_window_bp=(2 * SCORE_HALF + 1) * BIN, permissive_pad_bp=PERMISSIVE_PAD,
                            max_lag_bp=MAX_LAG_BINS * BIN)
            write_json(root / "thresholds.json", dict(**settings, **(thresholds or {})))
            rows = []
            with (root / "pooled_peaks.bed").open("w") as stream:
                for chrom, length in lengths.items():
                    original, short = geometry[chrom]
                    features = build_features(labs, cache, hard, permissive, chrom, length) if thresholds else None
                    for evidence in pool_chrom(original, short, envelope.get(chrom, []), length, features, thresholds):
                        name = f"pooled_peak_{len(rows) + 1:05d}"
                        rows.append(dict(chrom=chrom, peak=name, **evidence))
                        stream.write(f"{chrom}\t{evidence['start']}\t{evidence['end']}\t{name}\n")
        write_tsv(root / "pooled_peak_evidence.tsv", EVIDENCE_FIELDS, rows)
        summary = dict(status=settings["status"], peaks=len(rows), span_bp=sum(row["bp"] for row in rows),
                       **{field: sum(row[field] for row in rows) for field in EVIDENCE_FIELDS[5:]})
        if thresholds is None:
            summary["reason"] = "No grid setting had defined shifted support at or below the target on calibration chromosomes"
        write_json(root / "pooling_summary.json", summary)
    return summary
