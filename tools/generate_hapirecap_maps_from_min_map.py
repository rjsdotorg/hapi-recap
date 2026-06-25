#!/usr/bin/env python3
"""Generate HAPI-RECAP-compatible map files from a BIM and min_map.txt.

HAPI-RECAP expects map files with columns:
    chrom  marker_id  genetic_position_cm  physical_position

This script interpolates the cM value for each BIM marker using the non-zero
positions from a Visual Phaser-style min_map.txt file.

If only one source map is available, the same interpolated map can be written to
sex-averaged, male, and female outputs so HAPI-RECAP at least has meaningful
genetic distances instead of all-zero PLINK map values.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _pick_column(columns, aliases):
    normalized = {str(col).strip().lower(): col for col in columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def load_min_map(path: Path) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    raw = pd.read_csv(path, sep="\t", header=0)
    chrom_col = _pick_column(raw.columns, ["chromosome", "chrom", "chr"])
    pos_col = _pick_column(raw.columns, ["position", "pos"])
    cm_col = _pick_column(raw.columns, ["cm", "length (cm)"])
    if chrom_col is None or pos_col is None or cm_col is None:
        raise KeyError("min_map.txt is missing required columns: Chromosome, Position, cM")

    raw = raw[[chrom_col, pos_col, cm_col]].copy()
    raw.columns = ["chrom", "position", "cm"]
    raw["chrom"] = pd.to_numeric(raw["chrom"], errors="coerce")
    raw["position"] = pd.to_numeric(raw["position"], errors="coerce")
    raw["cm"] = pd.to_numeric(raw["cm"], errors="coerce")
    raw = raw.dropna()
    raw["chrom"] = raw["chrom"].astype(int)
    raw = raw[raw["chrom"].between(1, 23)]

    by_pos = {}
    by_cm = {}
    for chrom, dfc in raw.groupby("chrom", sort=False):
        dfc = dfc.sort_values("position")
        by_pos[str(chrom)] = dfc["position"].to_numpy(dtype=np.float64)
        by_cm[str(chrom)] = dfc["cm"].to_numpy(dtype=np.float64)
    return by_pos, by_cm


def build_interpolated_map(bim_path: Path, dmap_positions: dict[str, np.ndarray], dmap_cms: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    with bim_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom = str(parts[0]).strip().upper().replace("CHR", "")
            if chrom == "X":
                chrom = "23"
            if not chrom.isdigit():
                continue
            c = int(chrom)
            if not (1 <= c <= 23):
                continue

            marker_id = str(parts[1])
            try:
                bp = int(parts[3])
            except ValueError:
                continue

            map_pos = dmap_positions.get(chrom)
            map_cm = dmap_cms.get(chrom)
            if map_pos is None or map_cm is None or len(map_pos) == 0:
                continue

            cm = float(np.interp(float(bp), map_pos, map_cm))
            rows.append((chrom, marker_id, cm, bp))

    return pd.DataFrame(rows, columns=["chrom", "marker_id", "cm", "bp"])


def write_map(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fout:
        for chrom, marker_id, cm, bp in df.itertuples(index=False):
            fout.write(f"{chrom}\t{marker_id}\t{cm:.6f}\t{bp}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bim", required=True, help="Path to BIM file")
    parser.add_argument("--min-map", required=True, help="Path to non-zero min_map.txt")
    parser.add_argument("--sex-avg-out", required=True, help="Output path for sex-averaged map")
    parser.add_argument("--male-out", required=True, help="Output path for male map")
    parser.add_argument("--female-out", required=True, help="Output path for female map")
    args = parser.parse_args()

    dmap_positions, dmap_cms = load_min_map(Path(args.min_map))
    interpolated = build_interpolated_map(Path(args.bim), dmap_positions, dmap_cms)
    if interpolated.empty:
        raise ValueError("No BIM markers could be interpolated from min_map.txt")

    write_map(interpolated, Path(args.sex_avg_out))
    write_map(interpolated, Path(args.male_out))
    write_map(interpolated, Path(args.female_out))

    print(f"Wrote {len(interpolated)} markers to:")
    print(f"  {args.sex_avg_out}")
    print(f"  {args.male_out}")
    print(f"  {args.female_out}")


if __name__ == "__main__":
    main()