"""Choose signal-rescue thresholds with within-chromosome shifted controls."""

from itertools import product

import numpy as np

from .evidence import aligned_shape, flags, normalize_profiles
from .inputs import BIN


WEAK_GRID = (0.5, 0.75, 1.0, 1.25)
STRONG_GRID = (2.0, 2.25, 2.5, 2.75, 3.0, 3.5)
SHAPE_GRID = (0.55, 0.60, 0.65, 0.70, 0.75, 0.80)
SHIFTS_KB = ((211, -389), (-227, 457), (517, 875), (-571, -821), (383, -679), (-419, 733))
GRID = tuple(product(enumerate(WEAK_GRID), enumerate(STRONG_GRID), enumerate(SHAPE_GRID)))
METRIC_COLUMNS = (
    "split", "anchor_lab", "weak", "strong", "min_corr", "observed_tiles", "observed_eligible_tiles",
    "valid_shifts", "shifted_mean_tiles", "shifted_min_tiles", "shifted_max_tiles", "shifted_observed_ratio",
)


def empty_counts(shifts):
    shape = (len(WEAK_GRID), len(STRONG_GRID), len(SHAPE_GRID), 3, len(shifts) + 1, 2)
    return np.zeros(shape, dtype=np.int64)


def count_sustained(positions, mask):
    hits = positions[mask]
    if len(hits) < 2:
        return 0
    contiguous = np.diff(hits) == 1
    return np.count_nonzero(np.r_[contiguous, False] | np.r_[False, contiguous])


def calibrate_chrom(features, original, shifts_kb):
    scores, peers, profiles, called = features
    old = flags(original, scores.shape[1], center_only=False)
    # Compare the same anchored centers for all seven displacements without wrapping.
    margin = max(abs(k) for pair in shifts_kb for k in pair) * 1000 // BIN
    available = (called.sum(axis=0) == 1) & ~old
    available[:margin] = available[-margin:] = False
    totals = empty_counts(shifts_kb)
    for i, anchor_calls in enumerate(called):
        positions = np.flatnonzero(available & anchor_calls)
        if not len(positions):
            continue
        j, k = (peer for peer in range(len(called)) if peer != i)
        own = scores[i, positions]
        anchor = normalize_profiles(profiles[i], positions)
        for shift_id, (left, right) in enumerate(((0, 0), *shifts_kb)):
            pos_a = positions + left * 1000 // BIN
            pos_b = positions + right * 1000 // BIN
            eligible = ~called[j, pos_a] & ~called[k, pos_b]
            eligible_count = eligible.sum()
            a, b = peers[j, pos_a], peers[k, pos_b]
            corr_a = aligned_shape(anchor, profiles[j], pos_a)
            corr_b = aligned_shape(anchor, profiles[k], pos_b)
            for ic, corr in enumerate(SHAPE_GRID):
                eligible_a = np.where(corr_a >= corr, a, -np.inf)
                eligible_b = np.where(corr_b >= corr, b, -np.inf)
                low, high = np.minimum(eligible_a, eligible_b), np.maximum(eligible_a, eligible_b)
                for iw, weak in enumerate(WEAK_GRID):
                    eligible_own = eligible & (own >= weak)
                    weak_both = low >= weak
                    for js, strong in enumerate(STRONG_GRID):
                        accepted = eligible_own & (weak_both | (high >= strong))
                        totals[iw, js, ic, i, shift_id] = (
                            count_sustained(positions, accepted), eligible_count,
                        )
    return totals


def summarize(counts, labs, split):
    for (iw, weak), (js, strong), (ic, corr) in GRID:
        grid_cell = counts[iw, js, ic]
        grouped = [("all", grid_cell.sum(axis=0)), *zip((lab.name for lab in labs), grid_cell)]
        for name, by_shift in grouped:
            observed, eligible = map(int, by_shift[0])
            denominators = by_shift[1:, 1]
            valid = int(np.count_nonzero(denominators))
            mean = minimum = maximum = ratio = None
            if eligible and valid == len(denominators):
                scaled = by_shift[1:, 0] * (eligible / denominators)
                mean, minimum, maximum = float(scaled.mean()), float(scaled.min()), float(scaled.max())
                ratio = mean / observed if observed else None
            yield dict(split=split, anchor_lab=name, weak=weak, strong=strong, min_corr=corr,
                       observed_tiles=observed, observed_eligible_tiles=eligible, valid_shifts=valid,
                       shifted_mean_tiles=mean, shifted_min_tiles=minimum, shifted_max_tiles=maximum,
                       shifted_observed_ratio=ratio)


def select_thresholds(rows, target):
    fields = ("weak", "strong", "min_corr")
    candidates = [row for row in rows if row["split"] == "calibration" and row["anchor_lab"] == "all"
                  and row["shifted_observed_ratio"] is not None
                  and row["shifted_observed_ratio"] <= target]
    if not candidates:
        return None
    best = max(candidates, key=lambda r: (r["observed_tiles"], *(r[field] for field in fields)))
    thresholds = {field: best[field] for field in fields}
    holdout = next((row for row in rows if row["split"] == "holdout" and row["anchor_lab"] == "all"
                    and all(row[key] == value for key, value in thresholds.items())), None)
    return dict(**thresholds, calibration=best, holdout=holdout)
