"""Read the sample sheet, chromosome splits, and scored peak inputs."""

import csv
import gzip
import hashlib
import json
import math
import re
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import NamedTuple


BIN = 250
SPLITS = ("calibration", "holdout", "output_only")
SAMPLE_COLUMNS = ("lab", "sample", "role", "broadpeak", "permissive_broadpeak", "bam", "bigwig")
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


@dataclass(frozen=True)
class Sample:
    name: str
    bam: Path | None
    bigwig: Path | None
    broadpeak: Path | None
    permissive: Path | None

    @property
    def signal(self):
        return self.bam or self.bigwig


@dataclass(frozen=True)
class Lab:
    name: str
    targets: tuple[Sample, ...]
    controls: tuple[Sample, ...]


class Peak(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str
    log_p: str = ""


def read_table(path, required):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"{path} needs TSV columns {', '.join(required)}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"{path} has duplicate column names")
        for row in reader:
            if None in row:
                raise ValueError(f"{path} line {reader.line_num} has too many fields")
            yield {key: (value or "").strip() for key, value in row.items()}


def read_genome(path):
    lengths, splits = {}, {}
    for row in read_table(path, ("chrom", "length", "split")):
        chrom, length, split = row["chrom"], int(row["length"]), row["split"]
        if not chrom or chrom in lengths or any(char.isspace() for char in chrom):
            raise ValueError(f"Invalid or duplicate chromosome name {chrom!r} in {path}")
        if length <= 0 or split not in SPLITS:
            raise ValueError(f"{chrom} in {path} needs a positive length and a split from {SPLITS}")
        lengths[chrom] = length
        splits[chrom] = split
    if not lengths:
        raise ValueError(f"No chromosomes in {path}")
    return lengths, splits


def file_at(sheet, value):
    return (sheet.parent / Path(value).expanduser()).resolve() if value else None


def read_samples(path):
    names, groups = set(), {}
    for row in read_table(path, SAMPLE_COLUMNS[:3]):
        lab, name, role = (row[key] for key in SAMPLE_COLUMNS[:3])
        if not ID.fullmatch(lab) or not ID.fullmatch(name) or role not in ("target", "igg"):
            raise ValueError(f"Invalid lab, sample, or role for row {row} in {path}")
        if lab == "all" or name in names:
            raise ValueError(f"Lab 'all' is reserved and sample IDs must be unique ({name})")
        names.add(name)
        bam, bigwig, primary, permissive = (
            file_at(path, row.get(key, "")) for key in ("bam", "bigwig", "broadpeak", "permissive_broadpeak")
        )
        if role == "target" and primary is None:
            raise ValueError(f"{name} needs a primary MACS3 broadPeak")
        if role == "igg" and (primary or permissive):
            raise ValueError(f"{name} is an IgG control and needs only signal")
        groups.setdefault(lab, {"target": [], "igg": []})[role].append(
            Sample(name, bam, bigwig, primary, permissive)
        )
    if len(groups) != 3 or any(not group["target"] for group in groups.values()):
        raise ValueError("The three-lab pooling rule requires exactly three labs with at least one target each")
    return tuple(Lab(name, tuple(group["target"]), tuple(group["igg"])) for name, group in groups.items())


def read_peaks(path, lengths, scored=False, prefix=""):
    records = []
    order = {chrom: index for index, chrom in enumerate(lengths)}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        for no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4 or (scored and len(fields) != 9):
                format_name = "nine-column MACS3 broadPeak" if scored else "BED4 or broadPeak"
                raise ValueError(f"{path}:{no} needs {format_name} columns")
            chrom, start, end, name = fields[:4]
            if chrom not in lengths:
                continue
            start, end = int(start), int(end)
            if not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f"{path}:{no} interval is outside the declared chromosome length")
            log_p = fields[7] if scored else ""
            if scored and (not math.isfinite(float(log_p)) or float(log_p) < 0):
                raise ValueError(f"{path}:{no} requires a finite nonnegative -log10(p) in column 8")
            records.append(Peak(chrom, start, end, prefix + name, log_p))
    return sorted(records, key=lambda peak: (order[peak.chrom], peak.start, peak.end))


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def write_tsv(path, columns, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_bed(path, peaks):
    with path.open("w") as stream:
        for peak in peaks:
            stream.write("\t".join(map(str, peak[:4])) + "\n")


@contextmanager
def stage_directory(destination):
    if destination.exists():
        raise FileExistsError(f"{destination} already exists; choose a fresh output directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{destination.name}-", dir=destination.parent) as staging:
        temporary = Path(staging) / "output"
        temporary.mkdir()
        yield temporary
        temporary.rename(destination)
