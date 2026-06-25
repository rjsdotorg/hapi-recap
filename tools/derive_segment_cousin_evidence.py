#!/usr/bin/env python3
"""Derive segment-level sibling-cousin evidence directly from HAPI + IBD feathers.

This utility is automation-only and does not consume any workbook/visual phasing files.
It extracts HAPI parent segments for a chromosome and measures overlap (marker count)
between each segment and sibling-cousin IBD intervals from chrN.feather.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd


def extract_parent_segments(parents_key, inf_hapi_data):
    """Return parent segments as marker-index intervals per chromosome."""
    segments = defaultdict(list)

    prev_chrom = None
    most_recent_p_site = None
    codes = inf_hapi_data[parents_key]["codes"]

    for marker_idx, (inf_codes, cur_chrom) in enumerate(zip(codes, inf_hapi_data["chr"])):
        if cur_chrom != prev_chrom:
            if prev_chrom is not None:
                assert most_recent_p_site is not None
                first_marker = most_recent_p_site - inf_hapi_data["chrstr"][prev_chrom]
                last_marker = inf_hapi_data["chrstr"][cur_chrom] - 1 - inf_hapi_data["chrstr"][prev_chrom]
                if last_marker - first_marker >= 100:
                    segments[str(prev_chrom)].append((int(first_marker), int(last_marker)))
            most_recent_p_site = None
            prev_chrom = cur_chrom

        if (
            inf_codes is not None
            and (inf_codes[0] == "P" or inf_codes[0] == "PC" or inf_codes[0] == "PA")
        ):
            if most_recent_p_site is not None:
                first_marker = most_recent_p_site - inf_hapi_data["chrstr"][cur_chrom]
                last_marker = marker_idx - 1 - inf_hapi_data["chrstr"][cur_chrom]
                if last_marker - first_marker >= 100:
                    segments[str(cur_chrom)].append((int(first_marker), int(last_marker)))
            most_recent_p_site = marker_idx

    assert most_recent_p_site is not None
    first_marker = most_recent_p_site - inf_hapi_data["chrstr"][prev_chrom]
    last_marker = marker_idx - 1 - inf_hapi_data["chrstr"][prev_chrom]
    if last_marker - first_marker >= 100:
        segments[str(prev_chrom)].append((int(first_marker), int(last_marker)))

    return segments


def overlap_len(a_start, a_end, b_start, b_end):
    lo = max(int(a_start), int(b_start))
    hi = min(int(a_end), int(b_end))
    if hi < lo:
        return 0
    return hi - lo + 1


def main():
    parser = argparse.ArgumentParser(description="Derive HAPI segment-level sibling-cousin overlap evidence")
    parser.add_argument("--hapi-json", required=True, help="Path to HAPI all.json")
    parser.add_argument("--parents-key", required=True, help="Parent key in all.json (e.g., FAM1:100001-FAM1:100002)")
    parser.add_argument("--ibd-dir", required=True, help="Directory with chrN.feather files")
    parser.add_argument("--chromosome", required=True, type=int, help="Chromosome number to analyze (1-23)")
    parser.add_argument("--siblings", required=True, help="Comma-separated siblings (e.g., Ray,Tom,Diane)")
    parser.add_argument("--cousins", required=True, help="Comma-separated cousins (e.g., Wendy)")
    parser.add_argument("--out", default="", help="Optional CSV output path")
    args = parser.parse_args()

    siblings = [s.strip() for s in args.siblings.split(",") if s.strip()]
    cousins = [c.strip() for c in args.cousins.split(",") if c.strip()]
    chrom_key = str(args.chromosome)

    with open(args.hapi_json, "r", encoding="utf-8") as handle:
        inf_hapi_data = json.load(handle)

    if args.parents_key not in inf_hapi_data:
        raise KeyError(f"parents key not found in all.json: {args.parents_key}")

    all_segments = extract_parent_segments(args.parents_key, inf_hapi_data)
    chrom_segments = all_segments.get(chrom_key, [])
    if not chrom_segments:
        raise ValueError(f"No extracted parent segments for chromosome {chrom_key}")

    feather_path = Path(args.ibd_dir) / f"chr{args.chromosome}.feather"
    if not feather_path.exists():
        raise FileNotFoundError(f"Missing feather: {feather_path}")
    df = pd.read_feather(feather_path)

    records = []
    for seg_idx, (seg_start, seg_end) in enumerate(chrom_segments, start=1):
        seg_len = seg_end - seg_start + 1

        for sib in siblings:
            sib_rows = df[
                ((df["id1"].astype(str) == sib) & (df["id2"].astype(str).isin(cousins)))
                | ((df["id2"].astype(str) == sib) & (df["id1"].astype(str).isin(cousins)))
            ]

            total_overlap_markers = 0
            cousin_overlap = defaultdict(int)
            for _, row in sib_rows.iterrows():
                id1 = str(row["id1"])
                id2 = str(row["id2"])
                cousin = id2 if id1 == sib else id1
                ov = overlap_len(seg_start, seg_end, int(row["start"]), int(row["end"]))
                if ov > 0:
                    total_overlap_markers += ov
                    cousin_overlap[cousin] += ov

            records.append(
                {
                    "chromosome": int(args.chromosome),
                    "segment_index": seg_idx,
                    "segment_start_marker": seg_start,
                    "segment_end_marker": seg_end,
                    "segment_length_markers": seg_len,
                    "sibling": sib,
                    "total_overlap_markers": total_overlap_markers,
                    "overlap_fraction": (total_overlap_markers / seg_len) if seg_len > 0 else 0.0,
                    "cousin_overlap_breakdown": ";".join(f"{k}:{v}" for k, v in sorted(cousin_overlap.items())),
                }
            )

    out_df = pd.DataFrame.from_records(records)
    out_df = out_df.sort_values(["segment_index", "sibling"]).reset_index(drop=True)

    print(f"Chromosome {args.chromosome} extracted segments: {len(chrom_segments)}")
    print(out_df.to_string(index=False))

    if args.out:
        out_df.to_csv(args.out, index=False)
        print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
