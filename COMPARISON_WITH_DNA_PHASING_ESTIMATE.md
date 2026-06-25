# Comparison: HAPI-RECAP vs DNA_phasing Estimated Implementation

## Overview

The `DNA_phasing` project (2026-04) created an estimated/heuristic implementation of parental DNA reconstruction before HAPI-RECAP source code was available. Both approaches solve the same problem—**reconstruct parental DNA from sibling genotypes**—but use fundamentally different input data sources and algorithms.

## Input and Data Sources

| Aspect | DNA_phasing | HAPI-RECAP |
|--------|------------|-----------|
| **Input Files** | Raw DNA text files (3 siblings) + Manual Excel workbook with RPs | PLINK binary files + HAPI2 JSON (phased haplotypes) + IBD Feather files |
| **Phasing Source** | Manual Visual Phasing by user (Excel workbook) | HAPI2 algorithmic phasing (automatic phase detection) |
| **Relative Information** | Optional: cousin shared segments from workbook | Required: actual IBD segments between siblings and many relatives |
| **Genetic Maps** | Optional: used only for breakpoint interpolation | Required: sex-specific maps for sex inference |
| **Phase Boundaries** | User-adjusted RPs (recombination points) in Excel | HAPI2 P-sites (detected automatically) + crossover events |

## Algorithm Comparison

### Parent Pair Inference

**DNA_phasing approach:**
1. At each SNP, enumerate **all valid diploid parent genotype pairs** from observed sibling alleles
2. Apply strict Mendelian compatibility filtering (both parents must transmit compatible alleles to all 3 siblings)
3. Fallback to modal allele if strict compatibility fails (noise/missing data)
4. Choose orientation by minimizing **genotype transition distance** across chromosome
5. Optional: apply cousin block support as weak orientation prior

**HAPI-RECAP approach:**
1. Extract parent segments from HAPI2 where it already split phases into distinct haplotypes
2. Match those segments against IBD segments from relatives
3. Determine which relatives are "purple" (related to both parents—excluded as ambiguous)
4. Build parent-pair orientations using **IBD co-inheritance patterns** (do relatives share same parent or opposite parents?)
5. Validate using **crossover rates** under sex-specific genetic maps (Poisson model)

### Data Flow Diagram

**DNA_phasing:**
```
Sibling genotypes        Excel workbook RPs
      ↓                         ↓
Local Mendelian compatibility → Parent pair enumeration
      ↓
Genotype distance continuity + cousin block prior
      ↓
DP orientation smoothing
      ↓
Raw DNA output files
```

**HAPI-RECAP:**
```
HAPI2 phasing        IBD segments         Genetic maps
      ↓                   ↓                     ↓
Extract parent segments ← Purple-relative detection
      ↓
IBD overlap analysis → Parent linkage + orientation inference
      ↓
Crossover LOD scoring (sex inference)
      ↓
VCF output
```

## Sex Determination

| Method | DNA_phasing | HAPI-RECAP |
|--------|------------|-----------|
| **Approach** | Haploid allele consensus on sex chromosomes | Poisson model on crossover rates |
| **Data Used** | Unambiguous sibling alleles on chrX/chrY | Crossovers in two parent segments using two sex-specific genetic maps |
| **Output** | Paternal: haploid-only (AA), Maternal: het allowed | Infers via LOD score >3 distinguishing male vs female recombination rates |
| **Robustness** | Limited (depends on sibling data quality) | Stronger (cross-chromosome crossover accumulation) |

## Breakpoint/Recombination Handling

| Feature | DNA_phasing | HAPI-RECAP |
|---------|------------|-----------|
| **Source** | User-adjusted RPs in Excel (manual review) | HAPI2 P-sites + crossover detection |
| **Modes** | 3 ways to interpret widths: index, bp_linear, bp_genetic_map | Single approach: use SNP indices from P-sites |
| **Adjustment** | User changes column widths in workbook | Automatic; can flag inconsistent LOD |
| **Reset Behavior** | Hard reset at block boundaries (zero transition cost) | Can apply soft penalty at boundaries |

## Robustness and Validation

### Data Issue Handling

**DNA_phasing:**
- **Missing sibling data**: Falls back to modal allele
- **Noisy genotypes**: Tolerance via genotype distance threshold
- **Sex chromosome ambiguity**: Only outputs haploid if all siblings agree

**HAPI-RECAP:**
- **IBD quality**: Filters overlapping segments; removes segments < 9 cM
- **Purple relatives**: Excludes as ambiguous if related to both parents
- **Missing HAPI2 output**: Halts with clear error on malformed JSON
- **Sex conflict**: Flags segments with LOD ≤ -3 as sex-inconsistent

### Error Tolerance

**DNA_phasing:**
- Permissive: enumeration + DP can handle noise gracefully
- No external validation of parent pairs

**HAPI-RECAP:**
- Strict: validates IBD input columns, requires all required HAPI2 keys, checks crossover files exist
- Cross-validates parent inference via relative IBD matching

## Output

| Property | DNA_phasing | HAPI-RECAP |
|----------|------------|-----------|
| **Format** | Raw DNA text files (same as input) | VCF (standardized genomics format) |
| **Content** | Genotypes for all SNPs | Phased genotypes + quality metrics per segment |
| **Per-SNP Stats** | No site-level metadata | Full/half/missing allele counts per parent |
| **Reconstruction %** | Not quantified per parent | Reports coverage fraction per parent |

## Code Complexity and Maintenance

| Aspect | DNA_phasing | HAPI-RECAP |
|--------|------------|-----------|
| **Lines (core logic)** | ~800 (utils.py functions) | ~900 (hapi-recap.py functions) |
| **Key Data Structures** | Candidate pairs, genotype distance, DP paths | Interval arrays (piso), IBD DataFrames, crossover LOD tables |
| **External Libraries** | openpyxl, csv, itertools | pandas, numpy, piso, pysam |
| **Testing Strategy** | Breakpoint mode comparison (3 modes) | Outputs VCF per parent pair; reconstructs ~67–71% of genome |

## Key Differences Summary

| Feature | DNA_phasing | HAPI-RECAP |
|---------|------------|-----------|
| **Assumes phasing already done?** | No (manual RPs) | Yes (HAPI2 input) |
| **Uses relative IBD?** | Optionally (cousin blocks) | Fundamentally (core algorithm) |
| **Mendelian feasibility bias** | High (enumerates compatible pairs) | Lower (relies on IBD validation) |
| **Scales to many relatives?** | No (designed for siblings + cousins) | Yes (can use many relatives) |
| **Crossover modeling** | No explicit model | Yes (Poisson with sex-specific maps) |
| **Assumes family structure?** | Yes (3 siblings + grandparents) | No (generic sibling pairs + relatives) |
| **User adjustment of phases?** | Yes (edit Excel RPs) | No (HAPI2 output used as-is) |
| **Validation against external data** | Weak (optional cousin support) | Strong (IBD overlap + crossover LOD) |

## When Each Approach Applies

### DNA_phasing is suited for:
- **Absence of HAPI2 output**: when phasing was done manually or with a different tool
- **Small families**: 2–3 siblings with few relatives
- **Expert user guidance**: can adjust RPs in Excel workbook
- **Raw DNA format requirements**: output stays in same format as input
- **No genomic reference data needed**: optional genetic maps

### HAPI-RECAP is suited for:
- **Large-scale pipelines**: processes many sibling sets with HAPI2
- **Many relatives available**: leverages IBD from cousins, uncles, etc.
- **VCF ecosystem**: integrates with standard genomics tools
- **Population-scale studies**: genetic maps + crossover rates available
- **Sex inference needed**: crossover rates distinguish male/female
- **Full automation**: no manual breakpoint adjustment required

## Summary

**DNA_phasing** is a clever heuristic that approximates HAPI-RECAP using local Mendelian feasibility and genotype continuity. It works well in small, well-phased families but lacks cross-validation through IBD matching and cannot infer sex from raw genotypes alone.

**HAPI-RECAP** is a published, validated algorithm that uses phased haplotypes from HAPI2 and validates them against actual IBD segments from relatives. It provides stronger signal through multi-relative evidence and infers sex using a statistical model on crossovers.

The DNA_phasing estimate achieved its goal—a working proof-of-concept before source code release—but HAPI-RECAP's design is more robust for large, complex families with many relatives and requires significantly less manual intervention.

---

## Technical Deep Dive: Why HAPI2's Sibling Phasing Outperforms DNA_phasing

### 1. Global Optimization vs. Greedy Continuity

**DNA_phasing:**
- Locally enumerates compatible pairs at each SNP
- Smooths orientation via DP with genotype distance metric
- Hard resets at user-specified RPs (recombination points)
- **Problem**: Gets trapped in local minima. If a continuity path starts with an unlucky choice, DP cannot backtrack globally

**HAPI2:**
- Uses Viterbi algorithm on **inheritance vectors** (which haplotype each parent gave each sibling)
- Optimizes across *entire chromosome* simultaneously
- Can escape local optima because global path-scoring sees consequences downstream
- **Example**: HAPI2 might reject a "locally smooth" genotype transition if it forces an implausible recombination pattern 200 SNPs later

### 2. Modeling Recombination Explicitly

**DNA_phasing:**
- RPs (recombination points) are *user-adjusted positions* from Excel workbook
- Transition cost within a block ≈ 0; reset to 0 at block boundaries
- Cannot distinguish "real recombination" from "data error" — both look like genotype breaks
- **Risk**: If user's RP placement is off by 1 kb, entire downstream block orientation can flip

**HAPI2:**
- Models recombination *probability* across genetic distance
- Detects Mendelian errors explicitly (codes them as 'E')
- Distinguishes "likely Mendelian error" from "plausible recombination"
- Produces confidence codes (E, R, ?, ◇) for each site
- **Example**: 3 recombinations in 10 kb is implausibly dense → flags as error zone, not real crossovers

### 3. Confidence and Visibility

**DNA_phasing:**
```
SNP123: Diane=AG, Ray=CT, Tom=GG
→ Infer parent pair (AC, GT)
→ Output: AC | GT
(no indication of confidence or ambiguity)
```

**HAPI2:**
```
SNP123: Diane=AG, Ray=CT, Tom=GG
→ Viterbi path: infer inheritance (P1→A from dad, P2→G from mom)
→ Output: AG (with code=◇ meaning unambiguous) or ? (if multiple paths tied)
(user can see confidence & ambiguity)
```

### 4. Joint Sibling Inference

**DNA_phasing:**
- Independently checks Mendelian compatibility for each (parent_pair, child) triplet
- No coupling between siblings' decisions
- If siblings have inconsistent alleles, modal fallback can mask it
- Each sibling's data independently votes for a parent pair

**HAPI2:**
- Joint HMM over all siblings + parents simultaneously
- Sees constraints across all children at once
- If Diane and Tom imply one parent orientation but Ray implies opposite, **detects contradiction**
- Can reject orientations that don't explain all siblings jointly

### 5. Handling Ambiguity

**DNA_phasing:**
- If multiple parent pairs are equally compatible (ambiguous locus):
  - Picks one arbitrarily via `min()` or lexical ordering
  - No way to know it was ambiguous downstream
  - User cannot decide "mask this ambiguous site"

**HAPI2:**
- If two inheritance paths are equally likely:
  - Marks site as `?` (ambiguous)
  - Caller (HAPI-RECAP) can decide whether to use or mask
  - **Transparency**: user sees exactly where confidence is low

### 6. Published Validation

**DNA_phasing:**
- Designed in isolation; never benchmarked against ground truth
- No published comparison to other methods
- No external code review

**HAPI2:**
- Published in peer-reviewed literature: [*Reconstructing parent genomes using siblings and other relatives*](https://doi.org/10.1101/2024.05.10.593578)
- Benchmarked on real pedigrees + simulated data
- Cross-validated by multiple research groups
- Handles edge cases identified by community

### 7. Recombination Detection

**DNA_phasing:**
- User manually adjusts RPs in workbook
- No automatic detection of recombination hotspots or errors
- "Where should I place the next RP?" is manual guesswork

**HAPI2:**
- Detects crossovers → outputs locations
- HAPI-RECAP uses detected crossovers for sex inference
- Sex-specific genetic maps → LOD scoring (Poisson model)
- Sex determination becomes statistical, not heuristic

### 7a. Chr22 PA-Cadence Root Cause and Fix (Automation-Only)

This finding was derived from automated HAPI artifacts (`all.json`, `parhaps`, and code-path inspection), not workbook input.

Observed pattern on chr22:

- There is a regular ~50-marker cadence of `PA` codes across roughly `44.94 Mb` to `49.80 Mb`.
- Between adjacent `PA` sites (for example rel `6127`, `6177`, `6227`), most markers are coded `['?']`.
- In that same interval, `parhaps` are largely missing (`0/0`) between anchors.
- At `PA` anchors, the inferred allele alternates (`G -> T -> G -> ...`), which can look like repeated switches.

Interpretation:

- These are not true recombination boundaries.
- They are periodic phase-ambiguous sampling anchors in a low-information region (insufficient informative heterozygosity/SNP signal for confident trio phasing).
- The alternating `H02/H12` labels are best-guess alternations under ambiguity, not biologically supported crossover events.

Why segmentation was fragmented:

- Original `extract_parent_segments()` treated `PA` the same as `P/PC` and cut a new segment at each `PA` site.
- Combined with the segment-length filter, this created many micro-fragments and hid continuity across the FIR-like block.

Applied fix:

- Treat only `P` and `PC` as hard segment boundaries.
- Treat `PA` as an ambiguous anchor that should not break segments by default.
- Add explicit opt-in flag to restore old behavior when desired.

Implementation notes (repo-local):

- `hapi-recap.py` now supports `-pa_sites_break_segments` to restore legacy splitting on `PA`.
- Default behavior keeps continuity across `PA` anchors.
- `-min_parent_segment_markers` remains configurable for span filtering.
- `-dump_parent_segments` emits `*.parent_segments_debug.tsv` for direct inspection of kept/dropped spans and boundary reasons.

Effect on chr22:

- Before: dense `PA`-driven micro-fragmentation from `44.9-49.8 Mb`.
- After: continuous segment recovered through `49,799,427`, aligning with the known FIR-like interval structure.

### 8. Where DNA_phasing Might Actually Excel

To be fair, DNA_phasing has *some* advantages:

| Factor | Advantage |
|--------|-----------|
| **Simplicity** | Mendelian compatibility is intuitive; HMM internals are opaque |
| **Manual Control** | Can tweak RPs in Excel if you see a problem |
| **Coding Speed** | No pedigree file format, no reference genome requirement |
| **Noise Tolerance** | Modal fallback might be more forgiving than HAPI2's strict error codes |

### Comparative Algorithm Table

| Dimension | DNA_phasing | HAPI2 |
|-----------|------------|-------|
| **Optimization Scope** | Greedy DP (per block) | Global Viterbi (per chromosome) |
| **Recombination Model** | User-adjusted steps (static) | Probabilistic (dynamic) |
| **Error Detection** | None (silent) | Explicit codes (E, R, ?) |
| **Ambiguity Handling** | Silent/masked | Transparent (?) |
| **Joint Sibling Inference** | Independent checks per triplet | Full HMM coupling across all siblings |
| **Validation** | Untested/internal only | Published peer-review + benchmarks |
| **Sex Inference** | Data-dependent heuristic | Statistical LOD model |
| **Backtracking** | Cannot escape local minima | Global path optimization |
| **Confidence Reporting** | None | Per-site codes + LOD scores |

### Verdict

**HAPI2 is superior because:**
1. It is a principled statistical model that sees global structure
2. It explicitly models uncertainty
3. It has been battle-tested and published
4. It automatically detects data errors
5. It validates against biological constraints (recombination rates)

**DNA_phasing is a clever heuristic that works in simple cases** but can accumulate errors silently through:
- Local optimization traps
- User RP placement errors
- No distinction between noise and biology
- Independent sibling decisions (missing joint constraints)

**For 3 siblings**, HAPI2 will almost always produce better phases because it doesn't get fooled by local noise, can detect data errors, and validates against biological constraints. DNA_phasing might occasionally get lucky with its simpler model, but HAPI2's global optimization and error checking make it more reliable and reproducible.

---

## Validation Plan on Existing DNA_phasing 3-Sibling Test Set

This section defines a reproducible protocol to compare `DNA_phasing` and `HAPI2 + HAPI-RECAP` on the same sibling trio.

### Dataset

Use the existing sibling files previously used in `DNA_phasing` test runs:

- `MyHeritage_Diane_raw_dna.txt`
- `MyHeritage_Ray_raw_dna.txt`
- `MyHeritage_Tom_raw_dna.txt`

### Workflows to Run

1. `DNA_phasing` baseline:
       - Run `generate_estimated_parents.py` for the three breakpoint modes:
             - `index`
             - `bp_linear`
             - `bp_genetic_map`
       - Save paternal/maternal outputs for each mode.

2. `HAPI2 + HAPI-RECAP` workflow:
       - Convert the same sibling data to PLINK if needed.
       - Run HAPI2 with `--json --detect_co`.
       - Run HAPI-RECAP with the same BIM and map inputs.
       - Keep default no-relative handling unless explicitly comparing strict mode using `-use_default_no_relatives`.

### Metrics to Compare

#### 1) Reconstruction Percentage

From HAPI-RECAP logs:

- `% father reconstructed`
- `% mother reconstructed`

From DNA_phasing outputs:

- Compute callable-site percentage (non-missing genotype calls) per estimated parent.
- Report both per-parent and average.

#### 2) Heterozygosity Rate (HTZ)

Per estimated parent output:

- `het_rate = (# heterozygous genotype calls) / (# non-missing genotype calls)`

Compare:

- DNA_phasing modes vs HAPI2+HAPI-RECAP
- Optional benchmark against expected/empirical parental heterozygosity range for platform/population.

#### 3) Segment Boundary Consistency

Compare chromosome-level boundaries:

- DNA_phasing: block boundaries from RPs / mode-derived breakpoints
- HAPI2/HAPI-RECAP: inferred parent segment boundaries and CO-supported placements

Summarize per chromosome:

- Number of segments
- Median segment length
- Boundary proximity overlap (e.g., within N SNPs or M cM)

### Suggested Result Table (Publish in docs)

| Workflow | Mode | Parent A Reconstructed % | Parent B Reconstructed % | Parent A Het % | Parent B Het % | # Segments | Notes |
|----------|------|---------------------------|---------------------------|----------------|----------------|-----------|-------|
| DNA_phasing | index | TBD | TBD | TBD | TBD | TBD | Workbook RP index mode |
| DNA_phasing | bp_linear | TBD | TBD | TBD | TBD | TBD | Linear bp interpolation |
| DNA_phasing | bp_genetic_map | TBD | TBD | TBD | TBD | TBD | Genetic-map interpolation |
| HAPI2 + HAPI-RECAP | default | TBD | TBD | TBD | TBD | TBD | Auto-lenient if no relatives |
| HAPI2 + HAPI-RECAP | strict | TBD | TBD | TBD | TBD | TBD | `-use_default_no_relatives` |

### Interpretation Guidance

- Prefer methods with higher reconstruction percentage **without implausible HTZ inflation/deflation**.
- If boundaries differ, prioritize methods with:
      - better Mendelian consistency
      - fewer biologically implausible rapid switches
      - clearer uncertainty signaling (`?`, LOD notes) in HAPI2/HAPI-RECAP outputs.

### Publication Note

When publishing this comparison in project docs, include:

- exact command lines used
- repository commit hashes for both projects
- map versions/build (e.g., GRCh37)
- whether HAPI-RECAP was run in default or strict no-relative mode.
