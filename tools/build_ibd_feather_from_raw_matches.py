#!/usr/bin/env python
"""Build HAPI-RECAP ibd_chr*.feather directly from sibling/cousin raw DNA files.

This script ports the core Visual Phaser matching path (no workbook dependency):
- apply_conditions_vectorized()
- scan_genomes_optimized()

Workflow:
1) Load sibling and cousin raw DNA files.
2) For each sibling-cousin pair and chromosome, call HIR segments.
3) Convert segment BP bounds to marker-index bounds using BIM.
4) Write chr1..chr22 feather files in HAPI-RECAP schema.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_COLUMNS = ["id1", "id2", "chromosome", "start", "end", "start_cm", "end_cm"]


def _pick_column(columns, aliases):
    normalized = {str(col).strip().lower(): col for col in columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def _looks_like_vcf(file_path: Path) -> bool:
    if str(file_path).lower().endswith(".vcf"):
        return True
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for _ in range(20):
                line = handle.readline()
                if not line:
                    break
                stripped = line.strip()
                if stripped.startswith("##fileformat=VCF") or stripped.startswith("#CHROM"):
                    return True
    except OSError:
        return False
    return False


def _parse_vcf_file(file_path: Path) -> pd.DataFrame | None:
    header_columns = None
    separator = "\t"

    with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#CHROM"):
                if "\t" in line:
                    header_columns = line.lstrip("#").split("\t")
                    separator = "\t"
                else:
                    header_columns = line.lstrip("#").split()
                    separator = r"\s+"
                break

    if not header_columns:
        return None

    raw = pd.read_csv(
        file_path,
        sep=separator,
        comment="#",
        header=None,
        names=header_columns,
        dtype=str,
        low_memory=False,
        keep_default_na=False,
        engine="python" if separator != "\t" else None,
    )
    if raw.empty:
        return None

    chrom_col = _pick_column(raw.columns, ["chrom", "chromosome"])
    pos_col = _pick_column(raw.columns, ["pos", "position"])
    id_col = _pick_column(raw.columns, ["id", "rsid"])
    ref_col = _pick_column(raw.columns, ["ref"])
    alt_col = _pick_column(raw.columns, ["alt"])
    format_col = _pick_column(raw.columns, ["format"])

    if not all([chrom_col, pos_col, ref_col, alt_col]):
        return None

    df = pd.DataFrame(
        {
            "chromosome": raw[chrom_col].astype(str),
            "position": raw[pos_col].astype(str),
            "rsid": raw[id_col].astype(str) if id_col else "",
        }
    )

    ref_series = raw[ref_col].fillna("").astype(str)
    alt_series = raw[alt_col].fillna("").astype(str)
    df["allele1"] = ref_series
    df["allele2"] = alt_series.str.split(",").str[0]

    if format_col and len(header_columns) > 9:
        sample_col = header_columns[9]
        format_tokens = raw[format_col].fillna("").astype(str).str.split(":")
        gt_index = format_tokens.apply(lambda toks: toks.index("GT") if "GT" in toks else -1)
        sample_tokens = raw[sample_col].fillna("").astype(str).str.split(":")

        def _decode_gt(gt_idx, sample_vals, ref_val, alt_val):
            if gt_idx < 0 or gt_idx >= len(sample_vals):
                alt0 = alt_val.split(",")[0] if alt_val else ref_val
                return ref_val, alt0
            gt = sample_vals[gt_idx].replace("|", "/").strip()
            choices = [ref_val] + alt_val.split(",")
            parts = gt.split("/")

            def _pick(part):
                if not part or part == "." or not part.isdigit():
                    return ""
                idx = int(part)
                return choices[idx] if 0 <= idx < len(choices) else ""

            a1 = _pick(parts[0]) if len(parts) > 0 else ""
            a2 = _pick(parts[1]) if len(parts) > 1 else a1
            if not a1:
                a1 = ref_val
            if not a2:
                a2 = a1
            return a1, a2

        decoded = [
            _decode_gt(gt_idx, sample_vals, ref_val, alt_val)
            for gt_idx, sample_vals, ref_val, alt_val in zip(
                gt_index.tolist(),
                sample_tokens.tolist(),
                ref_series.tolist(),
                alt_series.tolist(),
            )
        ]
        df["allele1"] = [a for a, _ in decoded]
        df["allele2"] = [b for _, b in decoded]

    missing_id = df["rsid"].isin(["", ".", "nan", "None"])
    df.loc[missing_id, "rsid"] = (
        df.loc[missing_id, "chromosome"].astype(str)
        + ":"
        + df.loc[missing_id, "position"].astype(str)
    )
    return df


def _read_raw_table(file_path: Path) -> pd.DataFrame | None:
    if _looks_like_vcf(file_path):
        vcf_df = _parse_vcf_file(file_path)
        if vcf_df is not None:
            return vcf_df

    attempts = ["\t", ","]
    for sep in attempts:
        try:
            df = pd.read_csv(
                file_path,
                skip_blank_lines=True,
                comment="#",
                header=0,
                low_memory=False,
                dtype=str,
                keep_default_na=False,
                sep=sep,
            )
            if df is not None and len(df.columns) >= 4:
                return df
        except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError, OSError, ValueError):
            continue

    try:
        df = pd.read_csv(
            file_path,
            skip_blank_lines=True,
            comment="#",
            header=0,
            low_memory=False,
            dtype=str,
            keep_default_na=False,
            sep=None,
            engine="python",
        )
        if df is not None and len(df.columns) >= 4:
            return df
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError, OSError, ValueError):
        pass

    return None


def _clean_allele(series: pd.Series, no_call_val: str) -> pd.Series:
    no_call_token = str(no_call_val).strip().upper()
    cleaned = series.fillna("").astype(str).str.strip().str.upper()
    cleaned = cleaned.str.replace(r"[^A-Z0-9-]", "", regex=True)
    no_call_aliases = {"", "-", "--", "0", "00", "N", "NN", "NC", "NOCALL"}
    cleaned = cleaned.where(~cleaned.isin(no_call_aliases), no_call_token)
    cleaned = cleaned.where(cleaned.isin({"A", "T", "C", "G", no_call_token}), no_call_token)
    return cleaned


def load_individual_dna(ind: str, files_path: Path, no_call_val: str) -> pd.DataFrame:
    file_names = os.listdir(files_path)
    candidates = [name for name in file_names if f"{ind}_raw" in name]
    if not candidates:
        raise FileNotFoundError(f"No matching '*{ind}_raw*' file found in {files_path}")

    last_error = None
    for name in candidates:
        this_file = files_path / name
        try:
            raw = _read_raw_table(this_file)
            if raw is None or raw.empty:
                last_error = f"{name}: parser returned empty data"
                continue

            rsid_col = _pick_column(raw.columns, ["rsid", "rs#", "snp"])
            chrom_col = _pick_column(raw.columns, ["chromosome", "chrom", "chr"])
            pos_col = _pick_column(raw.columns, ["position", "pos"])
            allele1_col = _pick_column(raw.columns, ["allele1"])
            allele2_col = _pick_column(raw.columns, ["allele2"])
            genotype_col = _pick_column(raw.columns, ["result", "genotype", "alleles", "allele_pair"])

            if rsid_col is None or chrom_col is None or pos_col is None:
                cols = list(raw.columns)
                if len(cols) >= 4:
                    rsid_col, chrom_col, pos_col = cols[0], cols[1], cols[2]
                    if len(cols) >= 5:
                        allele1_col, allele2_col = cols[3], cols[4]
                    else:
                        genotype_col = cols[3]
                else:
                    last_error = f"{name}: missing required columns"
                    continue

            df = pd.DataFrame(
                {
                    "rsid": raw[rsid_col].astype(str),
                    "chromosome": raw[chrom_col].astype(str),
                    "position": raw[pos_col].astype(str),
                }
            )

            if allele1_col is not None and allele2_col is not None:
                df["allele1"] = raw[allele1_col]
                df["allele2"] = raw[allele2_col]
            elif genotype_col is not None:
                genotype = raw[genotype_col].fillna("").astype(str).str.strip().str.upper()
                genotype = genotype.str.replace(r"[^A-Z0-9-]", "", regex=True)
                df["allele1"] = genotype.str[0]
                df["allele2"] = genotype.str[1]
            else:
                last_error = f"{name}: no allele columns or genotype column"
                continue

            df["chromosome"] = df["chromosome"].str.strip().str.upper().str.replace("CHR", "", regex=False)
            df["chromosome"] = df["chromosome"].replace({"X": "23", "XY": "23", "MT": "M"})
            df = df[~df["chromosome"].isin(["Y", "M"])]
            df = df[df["chromosome"].str.isnumeric()]
            df["chromosome"] = df["chromosome"].astype(int)
            df = df[df["chromosome"].between(1, 23)]

            df["position"] = pd.to_numeric(df["position"], errors="coerce")
            df = df.dropna(subset=["position"])
            df["position"] = df["position"].astype(int)

            df["allele1"] = _clean_allele(df["allele1"], no_call_val)
            df["allele2"] = _clean_allele(df["allele2"], no_call_val)

            if df.empty:
                last_error = f"{name}: no usable autosomal rows"
                continue

            print(f"Loaded DNA file successfully: {name} ({ind})")
            return df.sort_values(by=["chromosome", "position"]).reset_index(drop=True)
        except Exception as exc:
            last_error = f"{name}: {exc}"

    raise RuntimeError(last_error or f"Could not parse input for {ind}")


def apply_conditions_vectorized(al1x, al2x, al1y, al2y, no_call_val):
    cond_nc = (al1x == no_call_val) | (al1y == no_call_val)
    cond_crimson = (al1x == al2x) & (al1y == al2y) & (al1x != al1y)
    cond_limegreen = ((al1x == al1y) & (al2x == al2y)) | ((al1x == al2y) & (al2x == al1y))

    res = np.full(al1x.shape, "yellow", dtype=object)
    res[cond_limegreen] = "limegreen"
    res[cond_crimson] = "crimson"
    res[cond_nc] = "limegreen"
    return res


def scan_genomes_optimized(dm, chrom, hir_cutoff, fir_cutoff, hir_snp_min, fir_snp_min, mm_dist, dmap_positions, dmap_cms):
    matches = dm["match"].values
    positions = dm["position"].values
    length = len(matches)

    dx, ds = [], []
    nmms = 0
    segflag = fflag = False
    stpos = pos = fstpos = fpos = nsnps = fsnps = mmpos = 0

    def get_dcm(start, end):
        stcm = np.interp(start, dmap_positions, dmap_cms)
        fincm = np.interp(end, dmap_positions, dmap_cms)
        return fincm - stcm

    for i in range(length):
        m, p = matches[i], positions[i]
        if not segflag:
            if m in ("yellow", "limegreen"):
                nsnps, segflag, stpos = 1, True, p
                if m == "limegreen":
                    fsnps, fstpos, fflag = 1, p, True
        else:
            if m in ("yellow", "limegreen"):
                nsnps += 1
                pos = p
                if fflag:
                    if m == "limegreen":
                        fsnps, fpos = fsnps + 1, p
                    else:
                        fflag = False
                        if fsnps > fir_snp_min:
                            dcm = get_dcm(fstpos, fpos)
                            if dcm > fir_cutoff:
                                ds.append({"Chr": chrom, "Start Mb": fstpos, "Finish Mb": fpos, "No. SNPs": fsnps, "Length (cM)": round(dcm, 1)})
                        fsnps = 0
                else:
                    if m == "limegreen":
                        fsnps, fstpos, fflag = 1, p, True
            else:
                if fflag:
                    if fsnps > fir_snp_min:
                        dcm = get_dcm(fstpos, fpos)
                        if dcm > fir_cutoff:
                            ds.append({"Chr": chrom, "Start Mb": fstpos, "Finish Mb": fpos, "No. SNPs": fsnps, "Length (cM)": round(dcm, 1)})
                    fflag, fsnps = False, 0

                nmms += 1
                if nmms == 1:
                    mmpos = p
                else:
                    if p - mmpos < mm_dist * 1000:
                        segflag, nmms = False, 0
                        if nsnps > hir_snp_min:
                            dcm = get_dcm(stpos, pos)
                            if dcm > hir_cutoff:
                                dx.append({"Chr": chrom, "Start Mb": stpos, "Finish Mb": pos, "No. SNPs": nsnps, "Length (cM)": round(dcm, 1)})
                        nsnps = 0
                    else:
                        nmms, mmpos = 1, p

    if segflag and nsnps > hir_snp_min:
        dcm = get_dcm(stpos, pos)
        if dcm > hir_cutoff:
            dx.append({"Chr": chrom, "Start Mb": stpos, "Finish Mb": pos, "No. SNPs": nsnps, "Length (cM)": round(dcm, 1)})
    if fflag and fsnps > fir_snp_min:
        dcm = get_dcm(fstpos, fpos)
        if dcm > fir_cutoff:
            ds.append({"Chr": chrom, "Start Mb": fstpos, "Finish Mb": fpos, "No. SNPs": fsnps, "Length (cM)": round(dcm, 1)})

    return pd.DataFrame(dx), pd.DataFrame(ds)


def repair_files_optimized(dm, fir_snp_min, mm_dist):
    matches, positions = dm["match"].values, dm["position"].values
    length = len(matches)
    firs = fir_snp_min // 2
    is_limegreen = matches == "limegreen"
    new_matches = matches.copy()

    for i in range(firs + 1, length - firs - 1):
        if matches[i] in ("crimson", "yellow"):
            if np.all(is_limegreen[i - firs:i]) and np.all(is_limegreen[i + 1:i + firs]):
                new_matches[i] = "limegreen"

    crimson_idx = np.where(new_matches == "crimson")[0]
    if len(crimson_idx) > 0:
        mm_dst = mm_dist * 1000
        for i in range(len(crimson_idx)):
            curr_pos = positions[crimson_idx[i]]
            isolated = True
            if i > 0 and curr_pos - positions[crimson_idx[i - 1]] <= mm_dst:
                isolated = False
            if i < len(crimson_idx) - 1 and positions[crimson_idx[i + 1]] - curr_pos <= mm_dst:
                isolated = False
            if isolated:
                new_matches[crimson_idx[i]] = "yellow"

    dm = dm.copy()
    dm["match"] = new_matches
    return dm


def _load_min_map(path: Path) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    raw = pd.read_csv(path, sep="\t", header=0)
    chrom_col = _pick_column(raw.columns, ["chromosome", "chrom", "chr"])
    pos_col = _pick_column(raw.columns, ["position", "pos"])
    cm_col = _pick_column(raw.columns, ["cm", "length (cm)"])
    if chrom_col is None or pos_col is None or cm_col is None:
        raise KeyError("min_map is missing required columns: Chromosome, Position, cM")

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
        by_pos[chrom] = dfc["position"].to_numpy(dtype=np.float64)
        by_cm[chrom] = dfc["cm"].to_numpy(dtype=np.float64)
    return by_pos, by_cm


def _load_bim_positions(bim_path: Path) -> dict[str, np.ndarray]:
    by_chrom: dict[str, list[int]] = defaultdict(list)
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
            if not (1 <= c <= 22):
                continue
            try:
                bp = int(parts[3])
            except ValueError:
                continue
            by_chrom[str(c)].append(bp)
    return {chrom: np.asarray(vals, dtype=np.int64) for chrom, vals in by_chrom.items()}


def _load_genetic_map(map_path: Path) -> dict[str, np.ndarray]:
    by_chrom: dict[str, list[float]] = defaultdict(list)
    with map_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom = str(parts[0]).strip().upper().replace("CHR", "")
            if chrom == "X":
                chrom = "23"
            if not chrom.isdigit():
                continue
            c = int(chrom)
            if not (1 <= c <= 22):
                continue
            try:
                cm = float(parts[2])
            except ValueError:
                continue
            by_chrom[str(c)].append(cm)
    return {chrom: np.asarray(vals, dtype=np.float64) for chrom, vals in by_chrom.items()}


def _to_marker_interval(chrom: str, start_bp: int, end_bp: int, bim_positions: dict[str, np.ndarray]) -> tuple[int, int] | None:
    positions = bim_positions.get(chrom)
    if positions is None or positions.size == 0:
        return None
    s_idx = int(np.searchsorted(positions, start_bp, side="left"))
    e_idx = int(np.searchsorted(positions, end_bp, side="right"))
    s_idx = min(max(s_idx, 0), int(positions.size))
    e_idx = min(max(e_idx, 0), int(positions.size))
    if e_idx <= s_idx:
        return None
    return s_idx, e_idx


def _cm_from_idx(chrom: str, idx: int, gmap: dict[str, np.ndarray]) -> float | None:
    cmap = gmap.get(chrom)
    if cmap is None or cmap.size == 0:
        return None
    clamped = min(max(idx, 0), int(cmap.size - 1))
    return float(cmap[clamped])


def _empty_chrom_df() -> pd.DataFrame:
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--files-path", required=True, help="Directory with *_raw* DNA files")
    parser.add_argument("--siblings", required=True, help="Comma-separated sibling IDs")
    parser.add_argument("--cousins", required=True, help="Comma-separated cousin IDs")
    parser.add_argument("--id1", required=True, help="Parent ID to assign as id1 for all cousin segments")

    parser.add_argument("--min-map", required=True, help="Visual Phaser min_map.txt path")
    parser.add_argument("--bim", required=True, help="BIM file used for HAPI2")
    parser.add_argument("--sex-avg-map", required=True, help="Sex-averaged map used by HAPI-RECAP")
    parser.add_argument("--out-dir", required=True, help="Output directory for chr1..chr22.feather")

    parser.add_argument("--no-call", default="--")
    parser.add_argument("--hir-cutoff", type=float, default=6.0)
    parser.add_argument("--fir-cutoff", type=float, default=3.0)
    parser.add_argument("--hir-snp-min", type=int, default=150)
    parser.add_argument("--fir-snp-min", type=int, default=80)
    parser.add_argument("--mm-dist", type=float, default=2.0)
    parser.add_argument("--repair", action="store_true", help="Enable Visual Phaser mismatch smoothing")
    args = parser.parse_args()

    files_path = Path(args.files_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    siblings = [s.strip() for s in args.siblings.split(",") if s.strip()]
    cousins = [c.strip() for c in args.cousins.split(",") if c.strip()]
    if len(siblings) < 1 or len(cousins) < 1:
        raise ValueError("Provide at least one sibling and one cousin")

    individuals = sorted(set(siblings + cousins))
    dna_cache: dict[str, pd.DataFrame] = {}
    for ind in individuals:
        dna_cache[ind] = load_individual_dna(ind, files_path, args.no_call)

    dmap_positions, dmap_cms = _load_min_map(Path(args.min_map))
    bim_positions = _load_bim_positions(Path(args.bim))
    gmap = _load_genetic_map(Path(args.sex_avg_map))

    raw_segments = []

    for chrom in range(1, 23):
        pos_arr = dmap_positions.get(chrom)
        cm_arr = dmap_cms.get(chrom)
        if pos_arr is None or cm_arr is None or len(pos_arr) == 0:
            continue

        for sib in siblings:
            for cousin in cousins:
                if sib == cousin:
                    continue

                dna1 = dna_cache[sib]
                dna2 = dna_cache[cousin]
                dm = pd.merge(
                    dna1[dna1["chromosome"] == chrom],
                    dna2[dna2["chromosome"] == chrom],
                    on=("rsid", "chromosome", "position"),
                    suffixes=("_1", "_2"),
                )
                if dm.empty:
                    continue

                dm = dm.sort_values("position").reset_index(drop=True)
                dm["match"] = apply_conditions_vectorized(
                    dm["allele1_1"].values,
                    dm["allele2_1"].values,
                    dm["allele1_2"].values,
                    dm["allele2_2"].values,
                    args.no_call,
                )

                if args.repair:
                    dm = repair_files_optimized(dm, args.fir_snp_min, args.mm_dist)

                dx, _ = scan_genomes_optimized(
                    dm,
                    chrom,
                    args.hir_cutoff,
                    args.fir_cutoff,
                    args.hir_snp_min,
                    args.fir_snp_min,
                    args.mm_dist,
                    pos_arr,
                    cm_arr,
                )
                if dx.empty:
                    continue

                for _, row in dx.iterrows():
                    start_bp = int(row["Start Mb"])
                    end_bp = int(row["Finish Mb"])
                    # Derive cM directly from min_map interpolation so cM lengths are
                    # meaningful even when external .map files have 0.0 genetic distances.
                    start_cm_interp = float(np.interp(start_bp, pos_arr, cm_arr))
                    end_cm_interp = float(np.interp(end_bp, pos_arr, cm_arr))
                    raw_segments.append(
                        {
                            "id1": str(sib),
                            "id2": str(cousin),
                            "chromosome": str(chrom),
                            "start_bp": start_bp,
                            "end_bp": end_bp,
                            "start_cm": start_cm_interp,
                            "end_cm": end_cm_interp,
                        }
                    )

    if not raw_segments:
        print("No segments found. Writing empty chr*.feather files.")

    seg_df = pd.DataFrame(raw_segments)

    for chrom in range(1, 23):
        chrom_str = str(chrom)
        out_path = out_dir / f"chr{chrom}.feather"

        if seg_df.empty:
            _empty_chrom_df().to_feather(out_path)
            continue

        chrom_rows = seg_df[seg_df["chromosome"] == chrom_str].copy()
        out_rows = []

        for _, row in chrom_rows.iterrows():
            idx_pair = _to_marker_interval(
                chrom_str,
                int(row["start_bp"]),
                int(row["end_bp"]),
                bim_positions,
            )
            if idx_pair is None:
                continue
            start_idx, end_idx = idx_pair

            start_cm = float(row.get("start_cm", np.nan))
            end_cm = float(row.get("end_cm", np.nan))
            if not np.isfinite(start_cm) or not np.isfinite(end_cm):
                start_cm = _cm_from_idx(chrom_str, start_idx, gmap)
                end_cm = _cm_from_idx(chrom_str, max(end_idx - 1, 0), gmap)
            if start_cm is None or end_cm is None:
                continue

            out_rows.append(
                {
                    "id1": str(row["id1"]),
                    "id2": str(row["id2"]),
                    "chromosome": chrom_str,
                    "start": int(start_idx),
                    "end": int(end_idx),
                    "start_cm": float(start_cm),
                    "end_cm": float(end_cm),
                }
            )

        if not out_rows:
            chrom_df = _empty_chrom_df()
        else:
            chrom_df = pd.DataFrame(out_rows, columns=OUTPUT_COLUMNS)
            chrom_df["id1"] = chrom_df["id1"].astype("string")
            chrom_df["id2"] = chrom_df["id2"].astype("string")
            chrom_df["chromosome"] = chrom_df["chromosome"].astype("string")
            chrom_df["start"] = chrom_df["start"].astype("int64")
            chrom_df["end"] = chrom_df["end"].astype("int64")
            chrom_df["start_cm"] = chrom_df["start_cm"].astype("float64")
            chrom_df["end_cm"] = chrom_df["end_cm"].astype("float64")
            chrom_df = chrom_df[chrom_df["end"] > chrom_df["start"]]
            chrom_df = chrom_df.sort_values(["id1", "id2", "start", "end"], kind="mergesort")

        chrom_df.to_feather(out_path)

    print(f"Wrote HAPI-RECAP feather files to: {out_dir}")


if __name__ == "__main__":
    main()
