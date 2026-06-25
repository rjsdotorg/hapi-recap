#!/usr/bin/env python
"""Build HAPI-RECAP ibd_chr*.feather files from an IBD segment table.

Expected output schema per chromosome file:
- id1 (string)
- id2 (string)
- chromosome (string)
- start (int64)    # marker index (0-based, half-open interval start)
- end (int64)      # marker index (0-based, half-open interval end)
- start_cm (float64)
- end_cm (float64)

This script supports two segment coordinate modes:
1) index: segment start/end are already marker indices
2) bp:    segment start/end are base-pair positions and are mapped to marker indices

For bp mode, provide both:
- --bim (same BIM used for HAPI2)
- --sex-avg-map (same sex-averaged map passed to hapi-recap.py)
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED_OUTPUT_COLUMNS = ["id1", "id2", "chromosome", "start", "end", "start_cm", "end_cm"]


def _normalize_sample_id(sample_id: str, id_style: str) -> str:
    """Normalize sample IDs to match HAPI keys or legacy dash format."""
    sid = str(sample_id)
    if id_style == "keep":
        return sid
    fam_dash = re.fullmatch(r"([^:\s]+)-(\d+)", sid)
    fam_colon = re.fullmatch(r"([^:\s]+):(\d+)", sid)
    if id_style == "hapi":
        if fam_colon:
            return sid
        if fam_dash:
            return f"{fam_dash.group(1)}:{fam_dash.group(2)}"
        return sid
    if id_style == "legacy":
        if fam_dash:
            return sid
        if fam_colon:
            return f"{fam_colon.group(1)}-{fam_colon.group(2)}"
        return sid
    return sid


def _normalize_chrom(value: object) -> str | None:
    text = str(value).strip().upper()
    if text.startswith("CHR"):
        text = text[3:]
    if text in {"X", "23"}:
        return "23"
    if text in {"Y", "24"}:
        return "24"
    if text in {"M", "MT", "25", "26"}:
        return None
    if text.isdigit():
        return str(int(text))
    return None


def _infer_sep(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".seg", ".txt"}:
        return "\t"
    if suffix == ".csv":
        return ","
    return None


def _read_segments(path: Path, sep: str | None) -> pd.DataFrame:
    if sep is None:
        return pd.read_csv(path, sep=None, engine="python")
    return pd.read_csv(path, sep=sep)


def _load_bim_positions(bim_path: Path) -> dict[str, np.ndarray]:
    by_chrom: dict[str, list[int]] = defaultdict(list)
    with bim_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom = _normalize_chrom(parts[0])
            if chrom is None or chrom not in {str(c) for c in range(1, 23)}:
                continue
            try:
                bp = int(parts[3])
            except ValueError:
                continue
            by_chrom[chrom].append(bp)
    return {chrom: np.asarray(vals, dtype=np.int64) for chrom, vals in by_chrom.items()}


def _load_genetic_map(map_path: Path) -> dict[str, np.ndarray]:
    by_chrom: dict[str, list[float]] = defaultdict(list)
    with map_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom = _normalize_chrom(parts[0])
            if chrom is None or chrom not in {str(c) for c in range(1, 23)}:
                continue
            try:
                cm = float(parts[2])
            except ValueError:
                continue
            by_chrom[chrom].append(cm)
    return {chrom: np.asarray(vals, dtype=np.float64) for chrom, vals in by_chrom.items()}


def _empty_chrom_df(chrom: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "id1": pd.Series(dtype="string"),
            "id2": pd.Series(dtype="string"),
            "chromosome": pd.Series(dtype="string"),
            "start": pd.Series(dtype="int64"),
            "end": pd.Series(dtype="int64"),
            "start_cm": pd.Series(dtype="float64"),
            "end_cm": pd.Series(dtype="float64"),
        }
    )


def _to_marker_intervals_from_bp(
    chrom: str,
    start_bp: pd.Series,
    end_bp: pd.Series,
    bim_positions: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    if chrom not in bim_positions:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    positions = bim_positions[chrom]
    if positions.size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    # Half-open interval in marker-index coordinates:
    # start_idx = first marker with pos >= start_bp
    # end_idx   = first marker with pos > end_bp
    start_idx = np.searchsorted(positions, start_bp.to_numpy(dtype=np.int64), side="left")
    end_idx = np.searchsorted(positions, end_bp.to_numpy(dtype=np.int64), side="right")

    start_idx = np.clip(start_idx, 0, positions.size)
    end_idx = np.clip(end_idx, 0, positions.size)
    return start_idx.astype(np.int64), end_idx.astype(np.int64)


def _cm_from_idx(chrom: str, idx: np.ndarray, gmap: dict[str, np.ndarray]) -> np.ndarray:
    cmap = gmap.get(chrom)
    if cmap is None or cmap.size == 0:
        return np.array([], dtype=np.float64)

    # idx is half-open boundary; for end boundary at len(cmap), clamp to last marker cm.
    clamped = np.clip(idx, 0, cmap.size - 1)
    return cmap[clamped]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--segments", required=True, help="Input IBD segment file (TSV/CSV)")
    parser.add_argument("--out-dir", required=True, help="Output directory for chr1..chr22 feather files")

    parser.add_argument("--sep", default=None, help="Input separator. Default: infer by extension")

    parser.add_argument("--id1-col", default="id1")
    parser.add_argument("--id2-col", default="id2")
    parser.add_argument(
        "--id-style",
        choices=["hapi", "legacy", "keep"],
        default="hapi",
        help=(
            "Sample ID output style for feather id1/id2. "
            "'hapi' uses FAM:ID (default), 'legacy' uses FAM-ID, 'keep' leaves values unchanged."
        ),
    )
    parser.add_argument("--chrom-col", default="chromosome")
    parser.add_argument("--start-col", default="start")
    parser.add_argument("--end-col", default="end")
    parser.add_argument("--start-cm-col", default="start_cm")
    parser.add_argument("--end-cm-col", default="end_cm")

    parser.add_argument(
        "--coords",
        choices=["index", "bp"],
        default="index",
        help="How to interpret --start-col/--end-col",
    )

    parser.add_argument("--bim", default=None, help="Required for --coords bp")
    parser.add_argument("--sex-avg-map", default=None, help="Required for --coords bp")

    parser.add_argument(
        "--id1-allow",
        default=None,
        help="Optional comma-separated allowlist for id1 (e.g. 100001,100002,FAM1:100001,FAM1:100002)",
    )

    args = parser.parse_args()

    seg_path = Path(args.segments)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sep = args.sep if args.sep is not None else _infer_sep(seg_path)
    segments = _read_segments(seg_path, sep=sep)

    needed = [args.id1_col, args.id2_col, args.chrom_col, args.start_col, args.end_col]
    if args.coords == "index":
        needed.extend([args.start_cm_col, args.end_cm_col])

    missing = [c for c in needed if c not in segments.columns]
    if missing:
        raise KeyError(f"Missing required input columns: {', '.join(missing)}")

    work = pd.DataFrame()
    work["id1"] = segments[args.id1_col].astype("string")
    work["id2"] = segments[args.id2_col].astype("string")
    work["id1"] = work["id1"].map(lambda x: _normalize_sample_id(x, args.id_style)).astype("string")
    work["id2"] = work["id2"].map(lambda x: _normalize_sample_id(x, args.id_style)).astype("string")
    work["chromosome"] = segments[args.chrom_col].map(_normalize_chrom).astype("string")

    # Keep autosomes only for this pipeline.
    work = work[work["chromosome"].isin([str(c) for c in range(1, 23)])].copy()

    if args.id1_allow:
        allow = {x.strip() for x in args.id1_allow.split(",") if x.strip()}
        work = work[work["id1"].isin(allow)].copy()

    if args.coords == "index":
        work["start"] = pd.to_numeric(segments.loc[work.index, args.start_col], errors="coerce").astype("Int64")
        work["end"] = pd.to_numeric(segments.loc[work.index, args.end_col], errors="coerce").astype("Int64")
        work["start_cm"] = pd.to_numeric(segments.loc[work.index, args.start_cm_col], errors="coerce")
        work["end_cm"] = pd.to_numeric(segments.loc[work.index, args.end_cm_col], errors="coerce")
    else:
        if not args.bim or not args.sex_avg_map:
            raise ValueError("--bim and --sex-avg-map are required for --coords bp")

        bim_positions = _load_bim_positions(Path(args.bim))
        gmap = _load_genetic_map(Path(args.sex_avg_map))

        start_bp = pd.to_numeric(segments.loc[work.index, args.start_col], errors="coerce")
        end_bp = pd.to_numeric(segments.loc[work.index, args.end_col], errors="coerce")

        starts = np.full(work.shape[0], -1, dtype=np.int64)
        ends = np.full(work.shape[0], -1, dtype=np.int64)

        for chrom, idxs in work.groupby("chromosome", sort=False).groups.items():
            chrom_start_bp = start_bp.loc[idxs]
            chrom_end_bp = end_bp.loc[idxs]
            s_idx, e_idx = _to_marker_intervals_from_bp(chrom, chrom_start_bp, chrom_end_bp, bim_positions)
            if s_idx.size == 0:
                continue
            starts[np.array(idxs, dtype=np.int64)] = s_idx
            ends[np.array(idxs, dtype=np.int64)] = e_idx

        work["start"] = pd.Series(starts, index=work.index, dtype="Int64")
        work["end"] = pd.Series(ends, index=work.index, dtype="Int64")

        start_cm_vals = np.full(work.shape[0], np.nan, dtype=np.float64)
        end_cm_vals = np.full(work.shape[0], np.nan, dtype=np.float64)

        for chrom, idxs in work.groupby("chromosome", sort=False).groups.items():
            idxs_arr = np.array(idxs, dtype=np.int64)
            s = work.loc[idxs, "start"].astype("int64").to_numpy()
            e = work.loc[idxs, "end"].astype("int64").to_numpy()

            s_cm = _cm_from_idx(chrom, s, gmap)
            # end is half-open boundary; map boundary to nearest existing marker cM.
            e_cm = _cm_from_idx(chrom, np.maximum(e - 1, 0), gmap)
            if s_cm.size == 0 or e_cm.size == 0:
                continue
            start_cm_vals[idxs_arr] = s_cm
            end_cm_vals[idxs_arr] = e_cm

        work["start_cm"] = pd.Series(start_cm_vals, index=work.index)
        work["end_cm"] = pd.Series(end_cm_vals, index=work.index)

    # Final filtering and typing.
    work = work.dropna(subset=["start", "end", "start_cm", "end_cm", "id1", "id2", "chromosome"]).copy()
    work["start"] = work["start"].astype("int64")
    work["end"] = work["end"].astype("int64")
    work["start_cm"] = work["start_cm"].astype("float64")
    work["end_cm"] = work["end_cm"].astype("float64")

    # Keep only non-empty positive intervals.
    work = work[work["end"] > work["start"]].copy()

    # Write chr1..chr22 files.
    for chrom in range(1, 23):
        chrom_str = str(chrom)
        out_path = out_dir / f"chr{chrom}.feather"
        chrom_df = work[work["chromosome"] == chrom_str][REQUIRED_OUTPUT_COLUMNS].copy()

        if chrom_df.empty:
            chrom_df = _empty_chrom_df(chrom)
        else:
            chrom_df["id1"] = chrom_df["id1"].astype("string")
            chrom_df["id2"] = chrom_df["id2"].astype("string")
            chrom_df["chromosome"] = chrom_df["chromosome"].astype("string")
            chrom_df = chrom_df.sort_values(["id1", "id2", "start", "end"], kind="mergesort")

        chrom_df.to_feather(out_path)

    print(f"Wrote feather files to: {out_dir}")


if __name__ == "__main__":
    main()
