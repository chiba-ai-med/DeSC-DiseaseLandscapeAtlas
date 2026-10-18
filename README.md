# DeSC-DiseaseLandscapeAtlas

Disease-family-specific Energy Landscape Analysis (ELA) on DeSC claims data.
51 disease families, full ICD-10 code space (7,000+ codes), exact PCA via
OnlinePCA.jl.

## Motivation

The prior whole-cohort ELA (`DeSC-OnlineXXX-Landscaper`) produced a landscape
whose low-energy basins reflected high-frequency lifestyle diseases. This
project rebuilds the analysis with **disease-family-specific index cohorts**,
yielding landscapes where **low-energy basins correspond to clinically severe
states** (Energy HR < 1 in the mortality/hospitalization Cox).

## Disease Families

51 diseases defined in `config/disease_families.tsv`. Preprocessing: `all_comorbid`
(all comorbidity codes, index disease excluded). PCA dimension: 7.

## Workflow

Each `.smk` file is independent. All share `workflow/common.smk`. Run in order:

```bash
COMMON="--configfile config/config.yaml --cores <n> --use-singularity --singularity-args '-B /home/koki/dev/DeSC-OnlineXXX-Landscaper/data:/home/koki/dev/DeSC-OnlineXXX-Landscaper/data'"

# 1. Define index cohorts
snakemake -s workflow/cohort.smk $COMMON

# 2. Build feature matrices (binary COO + zstd)
snakemake -s workflow/matrix.smk $COMMON

# 3. Run OnlinePCA (exact, out-of-core)
snakemake -s workflow/onlinepca.smk $COMMON

# 4. Run Landscaper (ELA: 128 patterns → basins → energy)
snakemake -s workflow/landscaper.smk $COMMON

# 5. Cox energy sanity check (death + hospitalization)
snakemake -s workflow/cox.smk $COMMON

# 6. Compile cross-disease summary
snakemake -s workflow/summary.smk $COMMON

# 7. Build per-patient scores table
snakemake -s workflow/patient.smk $COMMON

# 8. PCA pair density plots
snakemake -s workflow/pca_pairs.smk $COMMON

# 9. Per-disease analyses (KM, tagcloud, trajectory-to-death)
snakemake -s workflow/analysis.smk $COMMON

# 10. Figures (landscape atlas, KM comparison, code trajectory)
snakemake -s workflow/figures.smk $COMMON

# 11. Comparison with conventional risk scores (CCI, Elixhauser)
snakemake -s workflow/comparison.smk $COMMON
```

## Output Structure

```
output/{disease}/all_comorbid/
  cohort/            — index_cohort.tsv, feature_matrix.coo.zst, person_map.tsv
  pca/7/             — Scores.mm, Eigen_vectors.csv
  landscaper/7/      — Allstates, Basin, Energy
  cox/7/             — energy_sanity_check.tsv, energy_sanity_check_per_basin.tsv
  patient/           — patient_scores.tsv
  comparison/        — comorbidity_indices.tsv, risk_score_comparison.tsv
summary/
  atlas_summary.tsv            — PASS/FAIL across 51 diseases
  risk_comparison/             — CCI/Elixhauser comparison, energy correlates
plot/{disease}/                — KM, tagcloud, trajectory, PCA pairs
plot/{disease}/all_comorbid/7/ — landscape figures
plot/atlas_heatmap.png         — PASS/FAIL heatmap
```

## Container Images

- `ghcr.io/rikenbit/onlinepcajl:df65589` — OnlinePCA.jl
- `ghcr.io/chiba-ai-med/landscaper:main` — Landscaper
- `docker://koki/desc_onlinexxx_landscaper:20260602` — Julia + R analysis

## Requirements

- Snakemake >= 8.11.1
- Singularity/Apptainer
- Julia (OnlinePCA.jl, MatrixMarket)
- R (survival, data.table, ggplot2, survminer, ggwordcloud, scales)
