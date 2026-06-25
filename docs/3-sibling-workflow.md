# 3-Sibling Workflow Guide (HAPI2 + HAPI-RECAP)

This guide documents a practical workflow for reconstructing parents when you only have 3+ siblings and no informative external relatives.

## Scope

This workflow covers:

1. Preparing PLINK input for sibling-only analysis
2. Running HAPI2 on sibling-only pedigrees
3. Running HAPI-RECAP when no close relatives are available
4. Producing VCF output files for each family key in the HAPI2 JSON

## Inputs Required

- Sibling genotype data (PLINK `.bed/.bim/.fam`, or raw DNA that can be converted)
- Sex-averaged genetic map (`-sex_avg_map`)
- Male genetic map (`-male_map`)
- Female genetic map (`-female_map`)
- HAPI2 crossover directory (`-co_dir`)

Optional/empty-relatives case:

- IBD feather files (`-ibd_feather`) can be present but may have no close relatives passing HAPI-RECAP thresholds.

## Step 1: Prepare PLINK Files

If your data are already PLINK binary files (`data.bed`, `data.bim`, `data.fam`), continue to Step 2.

If starting from raw DNA files, convert to PLINK format with your preferred harmonization/conversion pipeline first. Ensure:

- Marker coordinates are consistent across siblings
- Allele coding is normalized
- `.fam` encodes the sibling-only pedigree with shared father and mother IDs

Sibling-only pedigree rule for HAPI2:

- For each sibling, use the same father ID and mother ID in `.fam`
- Parent genotype rows may be absent in sibling-only datasets

## Step 2: Run HAPI2 (Sibling-Only)

From the HAPI2 repository:

```bash
./hapi2 -p /path/to/data_prefix -o /path/to/hapi2_out --json --detect_co
```

Equivalent expanded form:

```bash
./hapi2 -g /path/to/data.bed -s /path/to/data.bim -i /path/to/data.bim -o /path/to/hapi2_out --json --detect_co
```

Expected outputs used by HAPI-RECAP:

- `/path/to/hapi2_out.json`
- Crossover files in output directory: `co-[parent_key].1` ... `co-[parent_key].22`

## Step 3: Prepare IBD Feather Inputs

HAPI-RECAP expects chromosome-specific feather files and currently reads 22 files via a `chr1 -> chrN` filename substitution pattern.

Example naming pattern:

- `ibd_chr1.feather`, `ibd_chr2.feather`, ..., `ibd_chr22.feather`

Required columns per file:

- `id1`, `id2`, `chromosome`, `start`, `end`, `start_cm`, `end_cm`

Sibling-only/no-relative case behavior:

- Files may exist but contain no segments that qualify as close relatives (`> 600 cM` after filtering)
- HAPI-RECAP now handles this gracefully and can place segments by crossover evidence

## Step 4: Run HAPI-RECAP

Standard run (auto-lenient when no close relatives are found):

```bash
python hapi-recap.py \
  -inf_hapi_json /path/to/hapi2_out.json \
  -ibd_feather /path/to/ibd_chr1.feather \
  -bim /path/to/data.bim \
  -sex_avg_map /path/to/map_sex_avg.txt \
  -male_map /path/to/map_male.txt \
  -female_map /path/to/map_female.txt \
  -co_dir /path/to/hapi2_out \
  -out /path/to/recap_vcf
```

Force original/default threshold behavior with no relatives:

```bash
python hapi-recap.py ... -use_default_no_relatives
```

What "auto-lenient" means:

- If no close relatives are found for a family, HAPI-RECAP lowers crossover-only segment assignment threshold from `|LOD| >= 3.0` to `|LOD| >= 2.5`
- If `-use_default_no_relatives` is set, it keeps `|LOD| >= 3.0`

## Step 5: Output Interpretation

Outputs are written per HAPI2 family key as:

- `/path/to/recap_vcf/[parent_key].vcf`

Each VCF includes:

- Reconstructed parent genotypes
- Phasing information
- Coverage summary in logs (fraction reconstructed for each parent)

In sibling-only/no-relative runs, expect:

- More dependence on crossover-only placement
- Potentially lower confidence than runs with multiple relatives

## Troubleshooting

### Missing IBD files

If any `ibd_chrN.feather` file is missing, HAPI-RECAP raises a clear `FileNotFoundError`.

### Missing required columns

If feather files are malformed, HAPI-RECAP raises a `KeyError` listing missing columns.

### Missing HAPI2 JSON keys

If JSON lacks required keys (`marker`, `physpos`, `chr`, `chrstr`), HAPI-RECAP raises a `KeyError`.

## Reproducibility Tips

- Keep a fixed map release and genome build for all runs
- Archive `.fam` pedigree definitions used for sibling-only setup
- Store exact command lines in a run log
- Version both HAPI2 and HAPI-RECAP commits for every result set

## Validation Against DNA_phasing Outputs

Use the metrics tool in this repository to compare:

- DNA_phasing estimated parent files (text format)
- HAPI-RECAP parent VCF output

Tool path:

- `tools/compare_3sibling_metrics.py`

Example run:

```bash
python tools/compare_3sibling_metrics.py \
  --dna-parent-a /path/to/Paternal_estimate_index.txt \
  --dna-parent-b /path/to/Maternal_estimate_index.txt \
  --recap-vcf /path/to/recap_vcf/123-456.vcf \
  --dna-label-a DNA_Paternal \
  --dna-label-b DNA_Maternal \
  --out-json /path/to/results/compare_metrics.json \
  --out-md /path/to/results/compare_metrics.md
```

The tool reports for each parent output:

- `total_sites`
- `callable_sites`
- `callable_percent`
- `het_sites`
- `het_percent_of_callable`
- `callable_segment_count`
- `callable_segment_median_bp_length`

Interpretation:

- Higher `callable_percent` indicates broader reconstruction coverage.
- `het_percent_of_callable` should be interpreted with biological plausibility in mind (avoid inflated/implausible heterozygosity).
- Segment count/median length helps compare fragmentation vs continuity.

## Current Environment Status (2026-06-24)

The following sibling-test assets were detected on the current machine:

- `C:\Users\rjs\AppData\Local\DNA_phasing\DNA_files\Estimated_Paternal_raw_dna_20260428_084146.txt`
- `C:\Users\rjs\AppData\Local\DNA_phasing\DNA_files\Estimated_Maternal_raw_dna_20260428_084146.txt`

HAPI2 executable was not detected under:

- `C:\Users\rjs\Documents\GitHub\williamslab\hapi2`

So end-to-end HAPI2+HAPI-RECAP rerun is pending HAPI2 build.

Build and run next:

```bash
cd C:/Users/rjs/Documents/GitHub/williamslab/hapi2/genetio
make
cd ..
make
```

Then run HAPI2 and HAPI-RECAP, and compare outputs with:

```bash
python tools/compare_3sibling_metrics.py \
  --dna-parent-a C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/Estimated_Paternal_raw_dna_20260428_084146.txt \
  --dna-parent-b C:/Users/rjs/AppData/Local/DNA_phasing/DNA_files/Estimated_Maternal_raw_dna_20260428_084146.txt \
  --recap-vcf C:/path/to/hapi-recap-output/[parent_key].vcf \
  --dna-label-a DNA_Paternal \
  --dna-label-b DNA_Maternal \
  --out-json C:/path/to/results/compare_metrics.json \
  --out-md C:/path/to/results/compare_metrics.md
```
