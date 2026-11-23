"""Run the versioned within-lab MSPC workflow."""

import json
import re
import subprocess

from .inputs import read_peaks, sha256, stage_directory, write_bed, write_json, write_tsv


MSPC_FLAGS = ["-r", "bio", "-w", "1e-2", "-s", "1e-5", "-g", "1e-5", "-a", "0.05",
              "-m", "lowest", "--excludeHeader"]


def version_of(executable):
    version = subprocess.check_output([executable, "--version"], text=True).strip()
    if not re.fullmatch(r"6\.1\.0(?:\.0)?", version):
        raise ValueError(f"MSPC v6.1.0 is required, found {version!r}")
    return version


def run_consensus(labs, lengths, destination, executable, threads):
    multi = any(len(lab.targets) > 1 for lab in labs)
    version = version_of(executable) if multi else None
    with stage_directory(destination) as root:
        for directory in ("mspc_inputs", "mspc", "peak_beds"):
            (root / directory).mkdir()
        receipts = {}
        summary = []
        for lab in labs:
            records = {sample.name: read_peaks(sample.broadpeak, lengths, scored=True)
                       for sample in lab.targets}
            primary = {sample.name: sha256(sample.broadpeak) for sample in lab.targets}
            if len(lab.targets) == 1:
                peaks, = records.values()
                method = "single_replicate"
            else:
                paths = []
                for sample, peaks in records.items():
                    path = root / "mspc_inputs" / f"{sample}.bed"
                    with path.open("w") as stream:
                        for peak in peaks:
                            score = "300" if float(peak.log_p) > 300 else peak.log_p
                            stream.write("\t".join(map(str, (*peak[:4], score))) + "\n")
                    paths.append(path)
                method = "MSPC_consensus"
                peaks = []
                # An empty replicate cannot contribute to an all-replicate overlap.
                if all(any(float(peak.log_p) >= 2 for peak in values) for values in records.values()):
                    folder = root / "mspc" / lab.name
                    command = [executable, "-i", *map(str, paths), *MSPC_FLAGS, "-c", str(len(paths)),
                               "-d", str(threads), "-o", str(folder)]
                    log = root / "mspc" / f"{lab.name}.log"
                    with log.open("w") as stream:
                        result = subprocess.run(command, cwd=root, stdout=stream, stderr=subprocess.STDOUT)
                    if result.returncode:
                        tail = log.read_text(errors="replace").splitlines()[-12:]
                        raise ValueError(f"MSPC exited with code {result.returncode} in {lab.name}\n"
                                         + "\n".join(tail))
                    native = folder / "ConsensusPeaks.bed"
                    if not native.is_file():
                        raise ValueError(f"MSPC did not write ConsensusPeaks.bed for {lab.name}")
                    peaks = read_peaks(native, lengths, prefix=f"{lab.name}_")
                    write_json(folder / "command.json", command)
            bed = root / "peak_beds" / f"{lab.name}.{method}.bed"
            write_bed(bed, peaks)
            entry = dict(lab=lab.name, method=method, peaks=len(peaks),
                         span_bp=sum(peak.end - peak.start for peak in peaks),
                         bed=str(bed.relative_to(root)), sha256=sha256(bed))
            summary.append(entry)
            receipts[lab.name] = dict(**entry, primary_sha256=primary)
        write_tsv(root / "peak_summary.tsv", tuple(summary[0]), summary)
        write_json(root / "step3.json", dict(mspc_version=version, mspc_flags=MSPC_FLAGS,
                                           lengths=lengths, labs=receipts))


def checked_consensus(labs, lengths, step3, primary_hashes):
    receipt = json.loads((step3 / "step3.json").read_text())
    if receipt["lengths"] != lengths or set(receipt["labs"]) != {lab.name for lab in labs}:
        raise ValueError("Step 3 chromosome lengths or lab names differ from the requested inputs")
    result = []
    for lab in labs:
        entry = receipt["labs"][lab.name]
        if entry["primary_sha256"] != {sample.name: primary_hashes[sample.name] for sample in lab.targets}:
            raise ValueError(f"{lab.name} primary peaks differ from Step 3; rerun Step 3 with these inputs")
        path = step3 / entry["bed"]
        if sha256(path) != entry["sha256"]:
            raise ValueError(f"{path} differs from the Step 3 result")
        result.append(read_peaks(path, lengths))
    return result
