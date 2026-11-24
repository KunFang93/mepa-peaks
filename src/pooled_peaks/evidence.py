"""Compute Poisson enrichment and lab-level signal shapes."""

import numpy as np
from scipy.special import xlogy

from .inputs import BIN


SCORE_HALF = 3
FLANK_INNER = 10
FLANK_OUTER = 42
SHAPE_HALF = 8
MAX_LAG_BINS = 2
PERMISSIVE_PAD = 500


def flags(intervals, n, center_only=True):
    diff = np.zeros(n + 1, dtype=np.int32)
    if intervals:
        rows = np.asarray(intervals, dtype=np.int64)
        if center_only:
            boundary = (rows + BIN - 1 - BIN // 2) // BIN
        else:
            boundary = np.column_stack((rows[:, 0] // BIN, (rows[:, 1] + BIN - 1) // BIN))
        np.add.at(diff, boundary[:, 0], 1)
        np.add.at(diff, boundary[:, 1], -1)
    return diff.cumsum()[:-1] > 0


def sustained(mask):
    adjacent = np.zeros_like(mask)
    adjacent[1:] |= mask[:-1]
    adjacent[:-1] |= mask[1:]
    return mask & adjacent


def window_sums(counts, half, with_flanks=False):
    prefix = np.r_[0, np.cumsum(counts, dtype=np.int64)]
    center = np.arange(len(counts))
    a = np.maximum(0, center - half)
    b = np.minimum(len(counts), center + half + 1)
    observed = prefix[b] - prefix[a]
    if not with_flanks:
        return observed
    left = np.maximum(0, center - FLANK_OUTER)
    inner_left = np.maximum(0, center - FLANK_INNER)
    inner_right = np.minimum(len(counts), center + FLANK_INNER + 1)
    right = np.minimum(len(counts), center + FLANK_OUTER + 1)
    size = (inner_left - left) + (right - inner_right)
    nearby = (prefix[inner_left] - prefix[left]) + (prefix[right] - prefix[inner_right])
    local = np.divide(nearby * (b - a), size, out=np.zeros(len(counts)), where=size > 0)
    return observed, local


def lab_signals(labs, counts, depth):
    scores, profiles = [], []
    for lab in labs:
        control_rate = None
        if lab.controls:
            control_counts = sum(window_sums(counts[s.name], SCORE_HALF) for s in lab.controls)
            control_rate = (control_counts + 1.0) / sum(depth[s.name] for s in lab.controls)
        scores_in_lab, shapes_in_lab = [], []
        for sample in lab.targets:
            name = sample.name
            observed, local = window_sums(counts[name], SCORE_HALF, with_flanks=True)
            mu = np.maximum(local, 0.5)
            if control_rate is not None:
                mu = np.maximum(mu, control_rate * depth[name])
            deviance = 2 * (xlogy(observed, observed / mu) - (observed - mu))
            scores_in_lab.append(np.sign(observed - mu) * np.sqrt(np.maximum(deviance, 0)))
            shapes_in_lab.append(window_sums(counts[name], 1) * (1e6 / depth[name]))
        scores.append(np.minimum.reduce(scores_in_lab))
        profiles.append(np.sqrt(np.mean(shapes_in_lab, axis=0)))
    return np.asarray(scores, dtype=np.float32), np.asarray(profiles, dtype=np.float32)


def normalize_profiles(signal, indices):
    windows = np.lib.stride_tricks.sliding_window_view(
        np.pad(signal, SHAPE_HALF, mode="edge"), 2 * SHAPE_HALF + 1
    )[indices]
    windows -= windows.mean(axis=1, keepdims=True)
    scale = np.linalg.norm(windows, axis=1, keepdims=True)
    np.divide(windows, scale, out=windows, where=scale > 0)
    return windows


def aligned_shape(anchor, profile, positions):
    best = np.full(len(positions), -np.inf, dtype=np.float32)
    for lag in range(-MAX_LAG_BINS, MAX_LAG_BINS + 1):
        shifted = np.clip(positions + lag, 0, len(profile) - 1)
        r = np.sum(anchor * normalize_profiles(profile, shifted), axis=1)
        np.maximum(best, r, out=best)
    return best


def candidate_anchors(scores, called, direct, min_labs):
    """Assign a partially called tile to its highest-scoring calling lab."""
    counts = called.sum(axis=0)
    anchors = np.full(len(counts), -1, dtype=np.intp)
    candidates = np.flatnonzero((counts > 0) & (counts < min_labs)
                                & ~flags(direct, len(counts), center_only=False))
    anchors[candidates] = np.where(called[:, candidates], scores[:, candidates], -np.inf).argmax(axis=0)
    return anchors, counts


def compare_peers(features, calling_lab, positions, anchor, offsets_kb):
    """Compare signal from noncalling labs, with optional within-chromosome shifts."""
    _, peers, profiles, called = features
    labs = np.array([i for i in range(len(called)) if i != calling_lab], dtype=np.intp)
    shifted = positions + np.asarray(offsets_kb)[:, None] * (1000 // BIN)
    noncallers = ~called[labs[:, None], positions]
    usable = noncallers & ~called[labs[:, None], shifted]
    # The already-called labs stay fixed. A shifted noncaller must still lack a call.
    eligible = np.all(~noncallers | usable, axis=0)
    values = np.where(usable, peers[labs[:, None], shifted], -np.inf)
    correlations = np.full_like(values, -np.inf)
    for row, lab in enumerate(labs):
        valid = eligible & np.isfinite(values[row])
        if valid.any():
            correlations[row, valid] = aligned_shape(anchor[valid], profiles[lab], shifted[row, valid])
    return values, correlations, eligible


def support_paths(weak_votes, needed, best_score, strong_cutoff):
    with_strong_peer = (weak_votes >= needed) & (best_score >= strong_cutoff)
    with_extra_weak_peer = weak_votes > needed
    return with_strong_peer, with_extra_weak_peer


def build_features(labs, cache, hard, permissive, chrom, length):
    scores, profiles = lab_signals(labs, cache.chrom(chrom), cache.depth)
    n = scores.shape[1]
    called = np.stack([flags(track.get(chrom, []), n) for track in hard])
    lab_support = []
    for lab in labs:
        samples = []
        for sample in lab.targets:
            padded = [(max(0, a - PERMISSIVE_PAD), min(length, b + PERMISSIVE_PAD))
                      for a, b in permissive[sample.name].get(chrom, [])]
            samples.append(flags(padded, n))
        lab_support.append(np.logical_and.reduce(samples))
    peers = np.where(np.stack(lab_support), scores, -np.inf)
    return scores, peers, profiles, called
