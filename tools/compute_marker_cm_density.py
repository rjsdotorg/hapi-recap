#!/usr/bin/env python3
import pandas as pd
import numpy as np

bim = r"C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/hapi2_trio/sibs3.bim"
mapf = r"C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/hapi2_trio/sex_avg.interp.map.txt"

rows = []
with open(bim, "r", encoding="utf-8", errors="ignore") as f:
    for ln in f:
        p = ln.split()
        if len(p) < 4:
            continue
        c = str(p[0]).strip().upper().replace("CHR", "")
        if c == "X":
            c = "23"
        if not c.isdigit():
            continue
        ci = int(c)
        if not (1 <= ci <= 23):
            continue
        rows.append((ci, int(p[3])))

bdf = pd.DataFrame(rows, columns=["chrom", "bp"])

mdf = pd.read_csv(mapf, sep="\t")
cols = {str(c).strip().lower(): c for c in mdf.columns}
chrom_col = cols.get("chromosome") or cols.get("chrom") or cols.get("chr")
pos_col = cols.get("position") or cols.get("pos") or cols.get("bp")
cm_col = cols.get("cm") or cols.get("length (cm)") or cols.get("genetic_position")
if chrom_col is None or pos_col is None or cm_col is None:
    # Fallback for PLINK-style .map-like text files without clear headers.
    # Expected order: chrom, snp_id, cm, bp
    raw = pd.read_csv(mapf, sep="\t", header=None)
    if raw.shape[1] >= 4:
        mdf = raw.iloc[:, [0, 3, 2]].copy()
        mdf.columns = ["chrom", "bp", "cm"]
    else:
        raise KeyError(f"Could not detect map columns in {mapf}; columns={list(mdf.columns)}")
else:
    mdf = mdf[[chrom_col, pos_col, cm_col]].copy()
    mdf.columns = ["chrom", "bp", "cm"]
mdf["chrom"] = pd.to_numeric(mdf["chrom"], errors="coerce")
mdf["bp"] = pd.to_numeric(mdf["bp"], errors="coerce")
mdf["cm"] = pd.to_numeric(mdf["cm"], errors="coerce")
mdf = mdf.dropna()

stats = []
for chrom, g in bdf.groupby("chrom", sort=True):
    mg = mdf[mdf["chrom"] == chrom].sort_values("bp")
    if mg.empty:
        continue
    bps = g["bp"].to_numpy(dtype=float)
    cms = np.interp(
        bps,
        mg["bp"].to_numpy(dtype=float),
        mg["cm"].to_numpy(dtype=float),
    )
    span_cm = float(cms.max() - cms.min()) if len(cms) > 1 else 0.0
    n = len(cms)
    markers_per_cm = (n / span_cm) if span_cm > 0 else float("nan")
    stats.append((chrom, n, span_cm, markers_per_cm, markers_per_cm * 7.0 if span_cm > 0 else float("nan")))

sdf = pd.DataFrame(stats, columns=["chrom", "markers", "span_cm", "markers_per_cm", "markers_for_7cm"])
autosomes = sdf[sdf["chrom"].between(1, 22)]

print("Autosomal markers_per_cm median:", round(float(autosomes["markers_per_cm"].median()), 2))
print("Autosomal markers_for_7cm median:", round(float(autosomes["markers_for_7cm"].median()), 1))

chr22 = sdf[sdf["chrom"] == 22]
print("\nChr22:")
print(chr22.to_string(index=False))

if not chr22.empty:
    mpc22 = float(chr22.iloc[0]["markers_per_cm"])
    print("\nChr22 conversions:")
    print("500 markers ~=", round(500.0 / mpc22, 2), "cM")
    print("7 cM ~=", round(7.0 * mpc22, 1), "markers")
