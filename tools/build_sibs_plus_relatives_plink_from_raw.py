#!/usr/bin/env python
"""Build PLINK input for siblings plus optional relative samples from raw DNA files.

This is a preparation utility for cousin-informed HAPI-RECAP runs:
- Siblings are placed in one nuclear family (shared father/mother IDs).
- Relatives (for example paternal cousins) are added as founders in a separate family.
- Marker set is anchored to sibling overlap; relative-missing markers are encoded as 0 0.

Input raw files are expected as tab-delimited with columns:
rsid, chromosome, position, result

Outputs:
- <out_prefix>.ped
- <out_prefix>.map
- <out_prefix>.bed/.bim/.fam (via plink2 --pedmap --make-bed)
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

VALID = {"A", "C", "G", "T"}


def norm_chrom(value: str) -> str | None:
    text = str(value).strip().upper()
    if not text:
        return None
    if text in {"X", "XY"}:
        return "23"
    if text == "Y":
        return "24"
    if text in {"MT", "M", "25", "26"}:
        return None
    if text.isdigit():
        return text
    return None


def norm_gt(gt: str) -> tuple[str, str] | None:
    if gt is None:
        return None
    alleles = [ch for ch in str(gt).strip().upper() if ch in VALID]
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    if len(alleles) < 2:
        return None
    return alleles[0], alleles[1]


def read_raw(path: Path) -> dict[tuple[str, int, str], tuple[str, str]]:
    calls: dict[tuple[str, int, str], tuple[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
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
            calls[(chrom, pos, rsid)] = gt
    return calls


def parse_sample_arg(values: list[str]) -> list[tuple[Path, str]]:
    samples: list[tuple[Path, str]] = []
    for raw in values:
        if ":" not in raw:
            raise ValueError(
                f"Invalid sample spec '{raw}'. Use --sibling path:id or --relative path:id"
            )
        path_text, sample_id = raw.split(":", 1)
        path = Path(path_text)
        sample_id = sample_id.strip()
        if not sample_id:
            raise ValueError(f"Missing sample id in sample spec '{raw}'")
        samples.append((path, sample_id))
    return samples


def chrom_key(chrom: str) -> int:
    try:
        return int(chrom)
    except ValueError:
        return 10**9


def build_map_loci(sibling_calls: list[dict[tuple[str, int, str], tuple[str, str]]]) -> list[tuple[str, int, str]]:
    loci = set(sibling_calls[0].keys())
    for calls in sibling_calls[1:]:
        loci &= set(calls.keys())
    return sorted(loci, key=lambda k: (chrom_key(k[0]), k[1], k[2]))


def write_ped_map(
    out_prefix: Path,
    siblings: list[tuple[str, dict[tuple[str, int, str], tuple[str, str]]]],
    relatives: list[tuple[str, dict[tuple[str, int, str], tuple[str, str]]]],
    father_id: str,
    mother_id: str,
) -> tuple[Path, Path]:
    sibling_calls = [calls for _, calls in siblings]
    loci = build_map_loci(sibling_calls)

    map_path = out_prefix.with_suffix(".map")
    ped_path = out_prefix.with_suffix(".ped")

    with map_path.open("w", encoding="utf-8", newline="") as map_file:
        for chrom, pos, rsid in loci:
            map_file.write(f"{chrom}\t{rsid}\t0\t{pos}\n")

    with ped_path.open("w", encoding="utf-8", newline="") as ped_file:
        for sample_id, calls in siblings:
            row = ["FAM1", sample_id, father_id, mother_id, "0", "-9"]
            for locus in loci:
                gt = calls.get(locus)
                row.extend([gt[0], gt[1]] if gt else ["0", "0"])
            ped_file.write(" ".join(row) + "\n")

        for sample_id, calls in relatives:
            row = ["FAMREL", sample_id, "0", "0", "0", "-9"]
            for locus in loci:
                gt = calls.get(locus)
                row.extend([gt[0], gt[1]] if gt else ["0", "0"])
            ped_file.write(" ".join(row) + "\n")

    return ped_path, map_path


def run_plink2_makebed(plink2_exe: Path, ped_path: Path, out_prefix: Path) -> None:
    cmd = [
        str(plink2_exe),
        "--pedmap",
        str(ped_path.with_suffix("")),
        "--make-bed",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sibling",
        action="append",
        required=True,
        help="Sibling sample as path:id (repeat for each sibling)",
    )
    parser.add_argument(
        "--relative",
        action="append",
        default=[],
        help="Relative sample as path:id (repeat for each relative, e.g. cousins)",
    )
    parser.add_argument("--father-id", default="100001")
    parser.add_argument("--mother-id", default="100002")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--plink2-exe", required=True)
    args = parser.parse_args()

    sibling_specs = parse_sample_arg(args.sibling)
    if len(sibling_specs) < 3:
        raise ValueError("At least 3 siblings are required")
    relative_specs = parse_sample_arg(args.relative)

    siblings: list[tuple[str, dict[tuple[str, int, str], tuple[str, str]]]] = []
    for path, sample_id in sibling_specs:
        siblings.append((sample_id, read_raw(path)))

    relatives: list[tuple[str, dict[tuple[str, int, str], tuple[str, str]]]] = []
    for path, sample_id in relative_specs:
        relatives.append((sample_id, read_raw(path)))

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ped_path, map_path = write_ped_map(
        out_prefix=out_prefix,
        siblings=siblings,
        relatives=relatives,
        father_id=args.father_id,
        mother_id=args.mother_id,
    )

    run_plink2_makebed(Path(args.plink2_exe), ped_path, out_prefix)

    print(f"Wrote: {map_path}")
    print(f"Wrote: {ped_path}")
    print(f"Wrote: {out_prefix}.bed/.bim/.fam")


if __name__ == "__main__":
    main()
