"""Convert BAM or count BigWig inputs to the same 250 bp fragment-count grid."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pyBigWig
import pysam

from .inputs import BIN, write_tsv


REQUIRED_FLAGS = 0x1 | 0x2 | 0x40
EXCLUDED_FLAGS = 0x4 | 0x8 | 0x100 | 0x200 | 0x400 | 0x800
MAX_EXACT_BIGWIG_COUNT = 2**24
CHUNK = 131072


def bin_slices(lengths):
    slices = {}
    offset = 0
    for chrom, length in lengths.items():
        end = offset + (length + BIN - 1) // BIN
        slices[chrom] = slice(offset, end)
        offset = end
    return slices, offset


def check_lengths(path, available, lengths):
    for chrom, length in lengths.items():
        if available.get(chrom) != length:
            raise ValueError(f"{path} has {chrom} length {available.get(chrom)}, expected {length}")


def check_signal(sample, lengths):
    if sample.signal is None:
        raise ValueError(f"{sample.name} needs a BAM or a 250 bp midpoint-count BigWig for Step 4")
    if not sample.signal.is_file():
        raise FileNotFoundError(f"{sample.name} selected input does not exist: {sample.signal}")
    if sample.bam is not None:
        with pysam.AlignmentFile(str(sample.bam), "rb") as stream:
            if not stream.is_bam:
                raise ValueError(f"{sample.bam} needs to be a BAM")
            check_lengths(sample.bam, dict(zip(stream.references, stream.lengths)), lengths)
    else:
        with pyBigWig.open(str(sample.bigwig)) as stream:
            if not stream.isBigWig():
                raise ValueError(f"{sample.bigwig} needs to be a BigWig")
            check_lengths(sample.bigwig, stream.chroms(), lengths)


def count_bam(path, lengths):
    slices, size = bin_slices(lengths)
    counts = np.zeros(size, dtype=np.uint32)
    seen = pairs = mapq_zero = 0
    pending = []
    with pysam.AlignmentFile(str(path), "rb", threads=2) as stream:
        if not stream.is_bam:
            raise ValueError(f"{path} needs to be a BAM")
        check_lengths(path, dict(zip(stream.references, stream.lengths)), lengths)
        offsets = [slices[chrom].start if chrom in slices else -1 for chrom in stream.references]
        ref_lengths = stream.lengths
        for read in stream.fetch(until_eof=True):
            seen += 1
            flag, tid = read.flag, read.reference_id
            if (flag & REQUIRED_FLAGS != REQUIRED_FLAGS or flag & EXCLUDED_FLAGS
                    or tid < 0 or tid != read.next_reference_id or offsets[tid] < 0):
                continue
            tlen = read.template_length
            if tlen == 0:
                continue
            midpoint = min(read.reference_start, read.next_reference_start) + abs(tlen) // 2
            if not 0 <= midpoint < ref_lengths[tid]:
                raise ValueError(f"{path} contains an out-of-bounds fragment midpoint")
            pending.append(offsets[tid] + midpoint // BIN)
            mapq_zero += read.mapping_quality == 0
            if len(pending) == CHUNK:
                np.add.at(counts, pending, 1)
                pairs += len(pending)
                pending.clear()
        if pending:
            np.add.at(counts, pending, 1)
            pairs += len(pending)
    if int(counts.sum()) != pairs:
        raise ValueError(f"{path} caused integer overflow in the midpoint counts")
    return counts, seen, pairs, mapq_zero


def count_bigwig(path, lengths):
    slices, size = bin_slices(lengths)
    counts = np.zeros(size, dtype=np.uint32)
    with pyBigWig.open(str(path)) as stream:
        if not stream.isBigWig():
            raise ValueError(f"{path} needs to be a BigWig")
        check_lengths(path, stream.chroms(), lengths)
        for chrom, length in lengths.items():
            chrom_counts = counts[slices[chrom]]
            previous_end = 0
            for start, end, value in stream.intervals(chrom) or ():
                if (start < previous_end or start % BIN or end != min(start + BIN, length)
                        or not float(value).is_integer() or not 0 <= value <= MAX_EXACT_BIGWIG_COUNT):
                    raise ValueError(
                        f"{path} at {chrom}:{start}-{end} needs one 250 bp bin per raw integer "
                        "midpoint count (except a clipped terminal bin); run export-midpoints"
                    )
                chrom_counts[start // BIN] = int(value)
                previous_end = end
    return counts, int(counts.sum())


def cache_sample(task):
    sample, lengths, directory = task
    if sample.bam is not None:
        counts, records, depth, mapq_zero = count_bam(sample.bam, lengths)
        kind = "bam"
    else:
        counts, depth = count_bigwig(sample.bigwig, lengths)
        records = mapq_zero = None
        kind = "midpoint_bigwig"
    if depth == 0:
        raise ValueError(f"{sample.name} has zero retained fragment midpoints on the selected chromosomes")
    np.save(directory / f"{sample.name}.npy", counts)
    return dict(sample=sample.name, source_type=kind, source=str(sample.signal),
                source_records=records, fragments=depth, mapq0_read1=mapq_zero)


class MidpointCounts:
    def __init__(self, samples, lengths, cache, report, workers):
        jobs = [(sample, lengths, cache) for sample in samples]
        if workers == 1:
            rows = list(map(cache_sample, jobs))
        else:
            with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
                rows = list(pool.map(cache_sample, jobs))
        self.depth = {row["sample"]: row["fragments"] for row in rows}
        self.slices, _ = bin_slices(lengths)
        self.counts = {sample.name: np.load(cache / f"{sample.name}.npy", mmap_mode="r")
                       for sample in samples}
        write_tsv(report, tuple(rows[0]), rows)

    def chrom(self, chrom):
        return {sample: track[self.slices[chrom]] for sample, track in self.counts.items()}


def export_bigwig(bam, lengths, destination):
    if destination.exists():
        raise FileExistsError(destination)
    counts, _, pairs, mapq_zero = count_bam(bam, lengths)
    if not pairs or counts.max() > MAX_EXACT_BIGWIG_COUNT:
        raise ValueError("BigWig export requires positive total depth and counts exactly representable as float32")
    slices, _ = bin_slices(lengths)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".bigwig-", dir=destination.parent) as directory:
        path = Path(directory) / destination.name
        with pyBigWig.open(str(path), "w") as stream:
            stream.addHeader(list(lengths.items()), maxZooms=0)
            for chrom, length in lengths.items():
                values = counts[slices[chrom]]
                nonzero = np.flatnonzero(values)
                for left in range(0, len(nonzero), CHUNK):
                    indices = nonzero[left:left + CHUNK]
                    starts = indices * BIN
                    stream.addEntries([chrom] * len(indices), starts.tolist(),
                                      ends=np.minimum(starts + BIN, length).tolist(),
                                      values=values[indices].astype(float).tolist())
        destination.hardlink_to(path)
    return pairs, mapq_zero
