#!/usr/bin/env python3
import pandas as pd

xls = pd.ExcelFile('C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/cousins/siblings_5.xlsx')

# Get Chr22
df = pd.read_excel(xls, 'Chr22', header=None)

print("="*100)
print("CHR22 DIANE-RAY SEGMENT DATA")
print("="*100)
print("\nRow 13 (headers):", df.iloc[13].tolist())
print("\nSegments (rows 14-18):")
for idx in range(14, 19):
    row = df.iloc[idx]
    if pd.notna(row[1]) and isinstance(row[1], (int, float)):
        chrom = int(row[1])
        start = int(row[2])
        end = int(row[3])
        snps = int(row[4])
        cm = row[5]
        cols6_11 = row[6:12].tolist()
        print(f"\n  Segment {idx-13}:")
        print(f"    Chr{chrom}: {start:,} - {end:,}")
        print(f"    Size: {cm} cM, {snps} SNPs")
        print(f"    Labels (cols 6-11): {cols6_11}")

print("\n" + "="*100)
print("GEDCOM SHEET (cousin names):")
print("="*100)
gedcom_df = pd.read_excel(xls, 'GEDCOMs')
print(gedcom_df.to_string())
