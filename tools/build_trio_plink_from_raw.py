#!/usr/bin/env python
"""Build a sibling-only trio PLINK dataset from 3 raw DNA files.

Inputs are expected in a 4-column format:
rsid, chromosome, position, result

Output files:
- <out_prefix>.ped
- <out_prefix>.map
- <out_prefix>.bed/.bim/.fam (via plink2)
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

VALID = {"A", "C", "G", "T"}


def norm_chrom(val: str) -> str | None:
    t = str(val).strip().upper()
    if not t:
        return None
    if t in {"X", "XY"}:
        return "23"
    if t == "Y":
        return "24"
    if t in {"MT", "M", "25", "26"}:
        return None
    if t.isdigit():
        return t
    return None


def norm_gt(gt: str) -> tuple[str, str] | None:
    if gt is None:
        return None
    a = [ch for ch in str(gt).strip().upper() if ch in VALID]
    if len(a) == 1:
        a = [a[0], a[0]]
    if len(a) < 2:
        return None
    return a[0], a[1]


def read_raw(path: Path) -> dict[tuple[str, int, str], tuple[str, str]]:
    d: dict[tuple[str, int, str], tuple[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for row in rdr:
            if not row or len(row) < 4:
                continue
            if row[0].strip().lower() == "rsid":
                continue
            rsid = row[0].strip()
            chrom = norm_chrom(row[1])
            if not rsid or chrom is None:
                continue
            try:
                pos = int(float(row[2]))
            except ValueError:
                continue
            gt = norm_gt(row[3])
            if gt is None:
                continue
            d[(chrom, pos, rsid)] = gt
    return d


def chrom_key(c: str) -> int:
    try:
        return int(c)
    except ValueError:
        return 10**9


def write_map_ped(
    out_prefix: Path,
    samples: list[str],
    per_sample: list[dict[tuple[str, int, str], tuple[str, str]]],
    father_id: str,
    mother_id: str,
) -> tuple[Path, Path]:
    common = set(per_sample[0].keys())
    for d in per_sample[1:]:
        common &= set(d.keys())

    loci = sorted(common, key=lambda k: (chrom_key(k[0]), k[1], k[2]))

    map_path = out_prefix.with_suffix(".map")
    ped_path = out_prefix.with_suffix(".ped")

    with map_path.open("w", encoding="utf-8", newline="") as mf:
        for chrom, pos, rsid in loci:
            mf.write(f"{chrom}\t{rsid}\t0\t{pos}\n")

    with ped_path.open("w", encoding="utf-8", newline="") as pf:
        for sid, d in zip(samples, per_sample):
            row = ["FAM1", sid, father_id, mother_id, "0", "-9"]
            for key in loci:
                a1, a2 = d[key]
                row.extend([a1, a2])
            pf.write(" ".join(row) + "\n")

    return ped_path, map_path


def run_plink2_makebed(plink2_exe: Path, ped: Path, out_prefix: Path) -> None:
    cmd = [
        str(plink2_exe),
        "--pedmap",
        str(ped.with_suffix("")),
        "--make-bed",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-a", required=True)
    ap.add_argument("--raw-b", required=True)
    ap.add_argument("--raw-c", required=True)
    ap.add_argument("--id-a", default="Diane")
    ap.add_argument("--id-b", default="Ray")
    ap.add_argument("--id-c", default="Tom")
    ap.add_argument("--father-id", default="100001")
    ap.add_argument("--mother-id", default="100002")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--plink2-exe", required=True)
    args = ap.parse_args()

    raw_paths = [Path(args.raw_a), Path(args.raw_b), Path(args.raw_c)]
    sample_ids = [args.id_a, args.id_b, args.id_c]
    per_sample = [read_raw(p) for p in raw_paths]

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ped_path, map_path = write_map_ped(
        out_prefix=out_prefix,
        samples=sample_ids,
        per_sample=per_sample,
        father_id=args.father_id,
        mother_id=args.mother_id,
    )

    run_plink2_makebed(Path(args.plink2_exe), ped_path, out_prefix)

    print(f"Wrote: {map_path}")
    print(f"Wrote: {ped_path}")
    print(f"Wrote: {out_prefix}.bed/.bim/.fam")


if __name__ == "__main__":
    main()
