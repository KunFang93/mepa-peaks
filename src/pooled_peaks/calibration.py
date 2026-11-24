"""Choose signal-rescue thresholds with within-chromosome shifted controls."""

from itertools import product

import numpy as np

from .evidence import candidate_anchors, compare_peers, normalize_profiles, support_paths
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


def resolve_shifts(lab_count, custom=None):
    peer_count = lab_count - 1
    if custom is None:
        # Preserve both historical displacements. Space extra peers at least 137 kb apart.
        return tuple(pair[:peer_count] + tuple(
            (-1)**trial * (1021 + 61 * trial + 137 * peer)
            for peer in range(peer_count - 2)
        ) for trial, pair in enumerate(SHIFTS_KB))
    shifts = tuple(map(tuple, custom))
    if (len(shifts) != len(SHIFTS_KB) or len(set(shifts)) != len(shifts)
            or any(len(row) != peer_count or 0 in row or len(set(row)) != peer_count for row in shifts)):
        raise ValueError(f"Provide six distinct --shift-kb sets with {peer_count} distinct nonzero offsets each")
    return shifts


def empty_counts(shifts, lab_count, min_labs):
    shape = (len(WEAK_GRID), len(STRONG_GRID), len(SHAPE_GRID),
             lab_count, len(shifts) + 1, min_labs - 1, 2)
    return np.zeros(shape, dtype=np.int64)


def count_sustained(positions, mask, strata, stratum_count):
    hits = np.flatnonzero(mask)
    if len(hits) < 2:
        return np.zeros(stratum_count, dtype=np.int64)
    contiguous = np.diff(positions[hits]) == 1
    sustained = np.r_[contiguous, False] | np.r_[False, contiguous]
    return np.bincount(strata[hits[sustained]], minlength=stratum_count)


def calibrate_chrom(features, direct, shifts_kb, min_labs):
    scores, _, profiles, called = features
    anchors, call_counts = candidate_anchors(scores, called, direct, min_labs)
    margin = max(abs(k) for shifts in shifts_kb for k in shifts) * 1000 // BIN
    anchors[:margin] = anchors[-margin:] = -1
    totals = empty_counts(shifts_kb, len(called), min_labs)
    displacements = ((0,) * (len(called) - 1), *shifts_kb)
    stratum_count = min_labs - 1
    for i in range(len(called)):
        positions = np.flatnonzero(anchors == i)
        if not len(positions):
            continue
        own = scores[i, positions]
        anchor = normalize_profiles(profiles[i], positions)
        needed = min_labs - call_counts[positions]
        strata = call_counts[positions] - 1
        for shift_id, shifts in enumerate(displacements):
            values, correlations, eligible = compare_peers(features, i, positions, anchor, shifts)
            eligible_count = np.bincount(strata[eligible], minlength=stratum_count)
            totals[:, :, :, i, shift_id, :, 1] = eligible_count
            for ic, corr in enumerate(SHAPE_GRID):
                matched = np.where(correlations >= corr, values, -np.inf)
                best_score = matched.max(axis=0)
                for iw, weak in enumerate(WEAK_GRID):
                    eligible_own = eligible & (own >= weak)
                    weak_votes = (matched >= weak).sum(axis=0)
                    for js, strong in enumerate(STRONG_GRID):
                        strong_path, weak_path = support_paths(weak_votes, needed, best_score, strong)
                        accepted = eligible_own & (strong_path | weak_path)
                        totals[iw, js, ic, i, shift_id, :, 0] = count_sustained(
                            positions, accepted, strata, stratum_count,
                        )
    return totals


def summarize(counts, labs, split):
    for (iw, weak), (js, strong), (ic, corr) in GRID:
        grid_cell = counts[iw, js, ic]
        grouped = [("all", grid_cell.sum(axis=0)), *zip((lab.name for lab in labs), grid_cell)]
        for name, by_shift in grouped:
            observed, eligible = map(int, by_shift[0].sum(axis=0))
            represented = by_shift[0, :, 1] > 0
            denominators = by_shift[1:, represented, 1]
            valid = int(np.count_nonzero(np.all(denominators > 0, axis=1))) if eligible else 0
            mean = minimum = maximum = ratio = None
            if eligible and valid == len(denominators):
                observed_by_count = by_shift[0, represented, 1]
                scaled = (by_shift[1:, represented, 0] * (observed_by_count / denominators)).sum(axis=1)
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
