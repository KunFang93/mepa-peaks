"""Command-line entry points for Step 3 and Step 4."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from .calibration import resolve_shifts
from .consensus import run_consensus, version_of
from .inputs import read_genome, read_peaks, read_samples, resolve_min_labs
from .midpoints import check_signal, export_bigwig
from .pooling import run_pooling


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def fraction(value):
    number = float(value)
    if not 0 <= number <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return number


def get_mspc(labs, program):
    if not any(len(lab.targets) > 1 for lab in labs):
        return None
    executable = shutil.which(program)
    if not executable:
        raise FileNotFoundError("Install MSPC v6.1.0 or pass --mspc /path/to/mspc")
    # subprocess runs MSPC from the stage directory.
    return str(Path(executable).resolve())


def check_inputs(labs, lengths, executable):
    if executable:
        print(f"MSPC {version_of(executable)}")
    for lab in labs:
        for target in lab.targets:
            peaks = read_peaks(target.broadpeak, lengths, scored=True)
            print(f"{lab.name}\t{target.name}\t{len(peaks)} primary peaks")
            if target.permissive is not None:
                read_peaks(target.permissive, lengths)
        for sample in (*lab.targets, *lab.controls):
            check_signal(sample, lengths)
            print(f"{lab.name}\t{sample.name}\t{'BAM' if sample.bam else 'midpoint BigWig'}\t{sample.signal}")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="mepa-peaks",
        description="MEPA (Multi-Experiment Peak Aggregation) pools broad peaks from two or more labs.",
    )
    stages = parser.add_subparsers(dest="stage", required=True)
    commands = {
        "check": "Check peak files, signal headers, and MSPC version",
        "step3": "Build within-lab consensus peak BEDs with MSPC v6.1.0",
        "step4": "Fit the shift-test cutoffs and export cross-lab peaks",
        "run": "Run Step 3 and Step 4",
        "export-midpoints": "Make a 250 bp raw fragment-midpoint count BigWig from BAM",
    }
    for stage, description in commands.items():
        command = stages.add_parser(stage, help=description, description=description)
        command.add_argument("--genome", required=True, type=Path, help="TSV with chrom, length and split")
        if stage == "export-midpoints":
            command.add_argument("--bam", required=True, type=Path)
            command.add_argument("--output", required=True, type=Path)
            continue
        command.add_argument("--samples", required=True, type=Path, help="TSV with labs, peaks and signal paths")
        if stage != "check":
            command.add_argument("--outdir", required=True, type=Path, help="New directory for this stage")
        if stage in ("check", "step3", "run"):
            command.add_argument("--mspc", default="mspc", help="Path to MSPC v6.1.0 (default: mspc on PATH)")
        if stage in ("step3", "run"):
            command.add_argument("--mspc-threads", type=positive, default=4)
        if stage == "step4":
            command.add_argument("--step3-dir", required=True, type=Path, help="Completed Step 3 directory")
        if stage in ("step4", "run"):
            command.add_argument("--workers", type=positive, default=2, help="Parallel BAM/BigWig readers (default: 2)")
            command.add_argument("--min-labs", type=positive,
                                 help="Minimum labs supporting a pooled peak (default: strict majority)")
            command.add_argument("--target-shift-ratio", type=fraction, default=0.10,
                                 help="Maximum shifted/observed calibration ratio (default: 0.10)")
            command.add_argument("--shift-kb", nargs="+", type=int, action="append", metavar="KB",
                                 help="One offset per other lab; repeat six times to replace the default shifts")
    return parser


def dispatch(args):
    lengths, splits = read_genome(args.genome.resolve())
    if args.stage == "export-midpoints":
        pairs, mapq_zero = export_bigwig(args.bam.resolve(), lengths, args.output.resolve())
        print(f"Wrote {args.output} with {pairs:,} midpoints, including {mapq_zero:,} MAPQ 0 read 1s")
        return
    labs = read_samples(args.samples.resolve())
    if args.stage in ("check", "step3", "run"):
        executable = get_mspc(labs, args.mspc)
    if args.stage == "check":
        check_inputs(labs, lengths, executable)
        return

    destination = args.outdir.resolve()
    if args.stage in ("step4", "run"):
        min_labs = resolve_min_labs(len(labs), args.min_labs)
        shifts = resolve_shifts(len(labs), args.shift_kb)
        if "calibration" not in splits.values():
            raise ValueError("Label at least one chromosome calibration in the genome TSV")
        for lab in labs:
            for sample in (*lab.targets, *lab.controls):
                check_signal(sample, lengths)

    if args.stage == "step3":
        run_consensus(labs, lengths, destination, executable, args.mspc_threads)
        print(f"Step 3 complete: {destination / 'peak_beds'}")
        return
    if args.stage == "run":
        step3 = destination / "step3_mspc"
        destination = destination / "step4_pool"
        if destination.exists():
            raise FileExistsError(destination)
        run_consensus(labs, lengths, step3, executable, args.mspc_threads)
        print(f"Step 3 complete: {step3 / 'peak_beds'}", flush=True)
    else:
        step3 = args.step3_dir.resolve()
    summary = run_pooling(
        labs, lengths, splits, step3, destination, args.workers, args.target_shift_ratio, shifts, min_labs
    )
    print(json.dumps(summary, indent=2))
    print(f"Step 4 complete: {destination / 'pooled_peaks.bed'}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        dispatch(args)
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"{parser.prog}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
