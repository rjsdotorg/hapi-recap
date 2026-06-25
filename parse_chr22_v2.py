#!/usr/bin/env python3
import pandas as pd

xls = pd.ExcelFile('C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/cousins/siblings_5.xlsx')
df = pd.read_excel(xls, 'Chr22', header=None)

print("="*120)
print("CHR22 DIANE-RAY SEGMENTS WITH ALL COLUMNS")
print("="*120)

# Row 13: section header
# Row 14: column headers for the FIR table
# Rows 15-18: segments

col_names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

print("\nRow 14 (Column Headers):")
for i, (col_name, val) in enumerate(zip(col_names, df.iloc[14].tolist())):
    print(f"  {col_name}: {val}")

print("\n\nSegments (rows 15-18) with ALL columns:")
for row_idx in range(15, 19):
    row = df.iloc[row_idx]
    if pd.notna(row[1]) and isinstance(row[1], (int, float)):
        chrom = int(row[1])
        start = int(row[2])
        end = int(row[3])
        snps = int(row[4])
        cm = row[5]
        print(f"\nSegment row {row_idx - 14}:")
        print(f"  Location: chr{chrom}:{start:,}-{end:,} ({cm} cM, {snps} SNPs)")
        print(f"  Full row:")
        for col_name, val in zip(col_names, row.tolist()):
            print(f"    {col_name}: {val}")

print("\n" + "="*120)
print("KEY INFO:")
print("  Wendy's family: Schumacher")
print("  Paternal grandparents: Koran/Schumacher")
print("  Columns H-L (indices 7-11) contain grandparent labels")
print("="*120)
