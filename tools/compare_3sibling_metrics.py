#!/usr/bin/env python
"""Compare DNA_phasing parent estimates with HAPI-RECAP VCF parent outputs.

This utility computes per-parent summary metrics for two workflows:
- DNA_phasing estimated parent raw DNA files (tab/comma delimited text)
- HAPI-RECAP reconstructed parent genotypes from one VCF

Metrics:
- total_sites
- callable_sites
- callable_percent
- het_sites
- het_percent_of_callable
- callable_segment_count
- callable_segment_median_bp_length

Callable segments are contiguous callable runs per chromosome in position order.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from dataclasses import dataclass
from pathlib import Path

try:
    from pysam import VariantFile
except ImportError:  # pragma: no cover - runtime dependency check
    VariantFile = None

VALID_ALLELES = {"A", "C", "G", "T", "I", "D"}
MISSING_TOKENS = {"", "--", "00", "NO CALL", "NOCALL", "NA", "N/A", "."}


@dataclass
class Site:
    chrom: str
    pos: int
    genotype: str


def normalize_chromosome(value: str) -> str:
    text = str(value).strip()
    up = text.upper()
    if up in {"X", "XY", "23"}:
        return "23"
    if up in {"Y", "24"}:
        return "24"
    return text


def normalize_genotype(text: str | None) -> str:
    if text is None:
        return "--"
    raw = str(text).strip().upper()
    if raw in MISSING_TOKENS:
        return "--"
    alleles = [ch for ch in raw if ch in VALID_ALLELES]
    if not alleles:
        return "--"
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    if len(alleles) > 2:
        alleles = alleles[:2]
    return "".join(sorted(alleles))


def is_callable(gt: str) -> bool:
    return gt != "--"


def is_het(gt: str) -> bool:
    return is_callable(gt) and len(gt) == 2 and gt[0] != gt[1]


def detect_delimiter(sample: str) -> str:
    if "\t" in sample:
        return "\t"
    return ","


def parse_raw_parent_file(path: Path) -> list[Site]:
    sites: list[Site] = []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delim = detect_delimiter(sample)
        reader = csv.reader(handle, delimiter=delim)
        for row in reader:
            if not row:
                continue
            first = str(row[0]).strip() if row[0] is not None else ""
            if first.startswith("#"):
                continue
            if len(row) < 4:
                continue

            chrom_val = str(row[1]).strip()
            pos_val = str(row[2]).strip()
            gt_val = str(row[3]).strip()

            if chrom_val.lower() == "chromosome" or pos_val.lower() == "position":
                continue

            try:
                pos = int(float(pos_val))
            except ValueError:
                continue

            sites.append(
                Site(
                    chrom=normalize_chromosome(chrom_val),
                    pos=pos,
                    genotype=normalize_genotype(gt_val),
                )
            )

    sites.sort(key=lambda s: (int(s.chrom) if s.chrom.isdigit() else 10**9, s.pos))
    return sites


def parse_recap_vcf(path: Path) -> dict[str, list[Site]]:
    if VariantFile is None:
        raise ImportError(
            "pysam is required to read HAPI-RECAP VCF outputs. Install with: pip install pysam"
        )

    sample_sites: dict[str, list[Site]] = {}
    with VariantFile(str(path), "r") as vcf:
        samples = list(vcf.header.samples)
        if len(samples) != 2:
            raise ValueError(f"Expected exactly 2 parent samples in VCF, found {len(samples)}")

        for sample in samples:
            sample_sites[sample] = []

        for rec in vcf:
            chrom = normalize_chromosome(rec.chrom)
            pos = int(rec.pos)
            for sample in samples:
                call = rec.samples[sample]
                alleles = call.alleles
                if alleles is None or len(alleles) != 2 or alleles[0] is None or alleles[1] is None:
                    gt = "--"
                else:
                    gt = normalize_genotype("".join(alleles))
                sample_sites[sample].append(Site(chrom=chrom, pos=pos, genotype=gt))

    for sample in sample_sites:
        sample_sites[sample].sort(key=lambda s: (int(s.chrom) if s.chrom.isdigit() else 10**9, s.pos))
    return sample_sites


def callable_segments(sites: list[Site]) -> list[tuple[str, int, int]]:
    segments: list[tuple[str, int, int]] = []
    active_chrom: str | None = None
    start_pos: int | None = None
    prev_pos: int | None = None

    for site in sites:
        if is_callable(site.genotype):
            if active_chrom is None:
                active_chrom = site.chrom
                start_pos = site.pos
                prev_pos = site.pos
            elif site.chrom != active_chrom:
                assert start_pos is not None and prev_pos is not None
                segments.append((active_chrom, int(start_pos), int(prev_pos)))
                active_chrom = site.chrom
                start_pos = site.pos
                prev_pos = site.pos
            else:
                prev_pos = site.pos
        else:
            if active_chrom is not None:
                assert start_pos is not None and prev_pos is not None
                segments.append((active_chrom, int(start_pos), int(prev_pos)))
                active_chrom = None
                start_pos = None
                prev_pos = None

    if active_chrom is not None:
        assert start_pos is not None and prev_pos is not None
        segments.append((active_chrom, int(start_pos), int(prev_pos)))

    return segments


def summarize_sites(sites: list[Site]) -> dict[str, float | int | None]:
    total_sites = len(sites)
    callable_count = sum(1 for s in sites if is_callable(s.genotype))
    het_count = sum(1 for s in sites if is_het(s.genotype))

    segs = callable_segments(sites)
    lengths = [end - start for _chrom, start, end in segs]
    median_len = statistics.median(lengths) if lengths else None

    callable_pct = (100.0 * callable_count / total_sites) if total_sites else 0.0
    het_pct = (100.0 * het_count / callable_count) if callable_count else 0.0

    return {
        "total_sites": total_sites,
        "callable_sites": callable_count,
        "callable_percent": round(callable_pct, 4),
        "het_sites": het_count,
        "het_percent_of_callable": round(het_pct, 4),
        "callable_segment_count": len(segs),
        "callable_segment_median_bp_length": median_len,
    }


def markdown_table(results: dict) -> str:
    lines = [
        "| Workflow | Parent | total_sites | callable_sites | callable_% | het_sites | het_%_of_callable | segment_count | median_segment_bp |",
        "|----------|--------|------------:|---------------:|-----------:|---------:|------------------:|--------------:|------------------:|",
    ]

    for workflow in ("dna_phasing", "hapi_recap"):
        per_parent = results[workflow]["parents"]
        for parent_label, stats in per_parent.items():
            lines.append(
                "| {wf} | {parent} | {total} | {callable_sites} | {callable_pct} | {het_sites} | {het_pct} | {seg_count} | {median_bp} |".format(
                    wf=workflow,
                    parent=parent_label,
                    total=stats["total_sites"],
                    callable_sites=stats["callable_sites"],
                    callable_pct=stats["callable_percent"],
                    het_sites=stats["het_sites"],
                    het_pct=stats["het_percent_of_callable"],
                    seg_count=stats["callable_segment_count"],
                    median_bp=stats["callable_segment_median_bp_length"],
                )
            )

    return "\n".join(lines)


def build_results(
    dna_parent_a: Path,
    dna_parent_b: Path,
    recap_vcf: Path,
    dna_labels: tuple[str, str],
) -> dict:
    dna_a_sites = parse_raw_parent_file(dna_parent_a)
    dna_b_sites = parse_raw_parent_file(dna_parent_b)

    recap_sites = parse_recap_vcf(recap_vcf)
    recap_samples = list(recap_sites.keys())

    results = {
        "dna_phasing": {
            "inputs": {
                "parent_a": str(dna_parent_a),
                "parent_b": str(dna_parent_b),
            },
            "parents": {
                dna_labels[0]: summarize_sites(dna_a_sites),
                dna_labels[1]: summarize_sites(dna_b_sites),
            },
        },
        "hapi_recap": {
            "inputs": {"vcf": str(recap_vcf)},
            "parents": {
                recap_samples[0]: summarize_sites(recap_sites[recap_samples[0]]),
                recap_samples[1]: summarize_sites(recap_sites[recap_samples[1]]),
            },
        },
    }
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare DNA_phasing and HAPI-RECAP parent metrics")
    parser.add_argument("--dna-parent-a", required=True, help="DNA_phasing parent raw DNA file A")
    parser.add_argument("--dna-parent-b", required=True, help="DNA_phasing parent raw DNA file B")
    parser.add_argument("--recap-vcf", required=True, help="HAPI-RECAP output VCF with 2 parent samples")
    parser.add_argument("--dna-label-a", default="DNA_parent_A", help="Display label for DNA parent A")
    parser.add_argument("--dna-label-b", default="DNA_parent_B", help="Display label for DNA parent B")
    parser.add_argument("--out-json", default="", help="Optional path to write JSON summary")
    parser.add_argument("--out-md", default="", help="Optional path to write markdown table")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    results = build_results(
        dna_parent_a=Path(args.dna_parent_a),
        dna_parent_b=Path(args.dna_parent_b),
        recap_vcf=Path(args.recap_vcf),
        dna_labels=(args.dna_label_a, args.dna_label_b),
    )

    print(json.dumps(results, indent=2))
    print()
    md = markdown_table(results)
    print(md)

    if args.out_json:
        Path(args.out_json).write_text(json.dumps(results, indent=2), encoding="utf-8")
    if args.out_md:
        Path(args.out_md).write_text(md + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
