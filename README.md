HAPI-RECAP
==========

HAPI-RECAP is a tool that combines the phasing output from [HAPI2](https://github.com/williamslab/hapi2)
with IBD (identity-by-descent) segments to reconstruct the DNA of both parents of a set of siblings.
More detailed information is available in the [HAPI2 documentation](https://github.com/williamslab/hapi2) and the [HAPI-RECAP preprint](https://doi.org/10.1101/2024.05.10.593578).

## Overview

HAPI-RECAP reconstructs parental genomes from sibling genotypes when one or both parents are missing. It works by:

1. Taking phased haplotypes from HAPI2 (which phases siblings against each other or against partially/fully available parents)
2. Using IBD segments between siblings and their relatives to determine which phased haplotype belongs to which parent
3. Identifying crossover events to infer parent sexes using sex-specific genetic maps
4. Outputting VCF files with reconstructed parental genotypes

This approach is particularly useful for families where parents are deceased or genotypes are unavailable, but where siblings and other genetic relatives have genotype data available.

For a sibling-only practical runbook (3+ siblings, no informative relatives), see:
`docs/3-sibling-workflow.md`.

## Requirements

### From HAPI2

HAPI-RECAP requires output from a HAPI2 run on your data. When running HAPI2, you must include:

* `--json` or `--json_par`: Outputs phased haplotypes and metadata in JSON format that HAPI-RECAP reads
* `--detect_co`: Outputs detected crossover events for each chromosome (used to infer parent sexes)

Example HAPI2 invocation:
```bash
./hapi2 -p input_data -o hapi2_output --json --detect_co
```

This generates:
- `hapi2_output.json` - Contains phased haplotypes and metadata
- `hapi2_output/co-[siblings].[chromosome]` - Crossover location files for each family and chromosome (1-22)

### Data Files Required

HAPI-RECAP requires the following files as input:

1. **HAPI2 JSON output** (`-inf_hapi_json`): The JSON file produced by HAPI2 containing phased haplotypes and SNP metadata
2. **IBD segments** (`-ibd_feather`): IBD calls in Apache Feather format, one file per chromosome (e.g., `ibd_chr1.feather`, `ibd_chr2.feather`, etc.)
   - Format: DataFrame with columns `id1`, `id2`, `chromosome`, `start`, `end`, `start_cm`, `end_cm`
   - `id1` and `id2` are individual IDs; `id1` should be the parent ID and `id2` the relative
   - `start`/`end`: SNP indices on the chromosome
   - `start_cm`/`end_cm`: Genetic positions in centiMorgans
3. **BIM file** (`-bim`): PLINK BIM file used to run HAPI2 (contains marker IDs and alleles)
4. **Genetic maps** (`-sex_avg_map`, `-male_map`, `-female_map`): Plain-text genetic map files with columns:
   - Chromosome ID
   - SNP ID
   - Genetic position (cM)
   - Physical position (optional)
5. **Crossover directory** (`-co_dir`): Directory containing HAPI2 crossover detection files

### Dependencies

HAPI-RECAP requires Python 3 and the following packages:
- numpy
- pandas
- pysam
- piso

Install with:
```bash
pip install numpy pandas pysam piso
```

## Usage

### Running HAPI-RECAP

```bash
python hapi-recap.py \
  -inf_hapi_json hapi2_output.json \
  -ibd_feather ibd_chr1.feather \
  -bim input_data.bim \
  -sex_avg_map genetic_map_sex_avg.txt \
  -male_map genetic_map_male.txt \
  -female_map genetic_map_female.txt \
  -co_dir hapi2_output \
  -out output_vcfs
```

### Command-line Arguments

| Argument | Description |
|----------|-------------|
| `-inf_hapi_json` | Path to JSON file output by HAPI2 (use `--json` option when running HAPI2) |
| `-ibd_feather` | Path pattern to IBD Feather files; use `ibd_chr1.feather` and script will find `ibd_chr2.feather`, etc. |
| `-bim` | PLINK BIM file used for HAPI2 (required for allele information) |
| `-sex_avg_map` | Sex-averaged genetic map file (3-column: chrom, marker_id, genetic_position_cm) |
| `-male_map` | Male-specific genetic map file (same format) |
| `-female_map` | Female-specific genetic map file (same format) |
| `-co_dir` | Directory containing HAPI2 crossover files (e.g., `co-[siblings].1`, `co-[siblings].2`, etc.) |
| `-out` | Output directory where VCF files will be written |
| `-use_default_no_relatives` | Optional flag: keep default crossover-only threshold when no close relatives are found; without this flag, HAPI-RECAP is slightly more lenient in sibling-only runs |
| `-no_relatives_lod_threshold` | Optional float (default `2.5`): crossover-only LOD threshold used when no close relatives are found |
| `-sibling_only_mode` | Optional flag for 3-sibling-only runs: if no close relatives are found, allows CO-only segment assignment even when sex-LOD is uninformative |
| `-emit_input_format_csv` | Optional flag: also write one reconstructed parent CSV per parent (`rsid,chromosome,position,result`) in the output directory |

## Output

### VCF Files

HAPI-RECAP outputs one VCF file per sibling family to the specified output directory:

- **Filename**: `[parent1_id]-[parent2_id].vcf`
- **Content**: Reconstructed phased genotypes for both parents across all chromosomes
- **Format**: Standard VCF format with:
  - Chromosome, position, alleles
  - Genotype phasing information for each parent
  - Missing data encoded as `./.`

The two parent columns indicate genotypes inferred for parent 1 and parent 2 respectively. The script determines parent sex (inferred father vs. mother) using sex-specific genetic maps and crossover information.

### Reconstruction Statistics

During the run, the script logs:

- **Reconstruction fraction**: Percentage of each parent's variants successfully reconstructed
- **Segment counts**: Number of segments inferred using IBD vs. crossovers only
- **Site counts**: For each parent, the number of:
  - Fully reconstructed sites (both alleles inferred)
  - Half-reconstructed sites (one allele inferred)
  - Missing sites

Example log output:
```
Reconstructed 67.3% of father's and 71.5% of mother's variants
from 18 segments using IBD and 2 segments using crossovers (CO):
    Father Full    Father Half    Mother Full    Mother Half
IBD    45.2%        12.1%          48.6%          14.2%
CO     2.1%         1.3%           3.1%           1.8%
```

## How It Works

1. **Phase Extraction**: Extracts regions where HAPI2 distinctly separated the two parental haplotypes
2. **IBD Matching**: Uses IBD segments between parents and their relatives to determine which phased chromosome belongs to which parent
3. **Purple Relative Identification**: Detects relatives with IBD to both parents (hints at misphasing or parental relationships)
4. **Parent Assignment**: Links parent segments to specific parents based on IBD overlap patterns with relatives
5. **Sex Inference**: Uses crossover rates in sex-specific genetic maps to infer parent sexes
6. **Genotype Reconstruction**: Outputs reconstructed parental genotypes in VCF format

## Example Workflow

```bash
# Step 1: Prepare data and run HAPI2
cd /path/to/hapi2
./hapi2 -p ../my_data/family_data -o ../my_data/hapi2_results --json --detect_co

# Step 2: Have IBD segments computed (using external method like IBIS, KING, etc.)
# Output should be in Feather format: ibd_chr1.feather, ibd_chr2.feather, ...

# Step 3: Obtain genetic maps (e.g., from HapMap or other sources)

# Step 4: Run HAPI-RECAP
cd /path/to/hapi-recap
python hapi-recap.py \
  -inf_hapi_json ../my_data/hapi2_results.json \
  -ibd_feather ../my_data/ibd_chr1.feather \
  -bim ../my_data/family_data.bim \
  -sex_avg_map ../genetic_maps/map_sexAvg.txt \
  -male_map ../genetic_maps/map_male.txt \
  -female_map ../genetic_maps/map_female.txt \
  -co_dir ../my_data/hapi2_results \
  -out ../my_data/reconstructed_parents
```

## Citation

If you use HAPI-RECAP, please cite:

Qiao Y, Jewett EM, McManus KF, Freeman WA, Curran JE, Williams-Blangero S, Blangero J, The 23andMe Research Team, Williams AL. Reconstructing parent genomes using siblings and other relatives. bioRxiv. https://doi.org/10.1101/2024.05.10.593578
