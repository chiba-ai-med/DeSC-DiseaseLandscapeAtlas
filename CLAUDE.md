# CLAUDE.md

## Project Overview

Disease-family-specific Energy Landscape Analysis on DeSC claims data. Unlike
the prior project (`../DeSC-OnlineXXX-Landscaper`) which ran a single
whole-cohort landscape, this project creates separate landscapes per disease
family (dementia, diabetes, heart failure/cancer) and explores preprocessing
strategies to obtain landscapes where low-energy basins reflect disease
progression / severity.

## Key Design Decisions

- **Index cohort design**: Each disease family defines an index cohort (patients
  with the index disease). The index disease codes themselves are excluded from
  the feature matrix.
- **Energy Cox sanity check**: Every setting is evaluated by whether
  `HR(Energy) < 1` (low energy = high risk). This is a go/no-go gate, not a
  final result.
- **Preprocessing exploration**: Multiple strategies (all comorbidities,
  drop-top-N frequent codes, complication-only codes) are compared
  systematically.

## Workflow Execution

Each `.smk` file is independent. All share `workflow/common.smk` which reads
the disease list from `config/disease_families.tsv` (51 diseases). Run in order:

```bash
snakemake -s workflow/cohort.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/matrix.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/onlinepca.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/landscaper.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/cox.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/summary.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/patient.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/pca_pairs.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/analysis.smk --configfile config/config.yaml --cores <n>
snakemake -s workflow/figures.smk --configfile config/config.yaml --cores <n>
```

## Wildcards

Fixed by `config/config.yaml` and `workflow/common.smk`:

- `{disease}`: 51 families from `config/disease_families.tsv`
  (hf_cancer is a legacy combined cohort, excluded from the standard pipeline)
- `{preproc}`: `all_comorbid` (fixed)
- `{dim}`: `7` (fixed)

## Architecture

- **Julia**: Sparse matrix I/O, OnlinePCA, Landscaper input preparation
- **R**: Cohort definition, Cox analysis, visualization
- **Snakemake**: Workflow orchestration with Singularity containers
- **Config**: `config/config.yaml` (paths, parameters),
  `config/disease_families.tsv` (ICD-10 definitions)

## Source Files

- `src/define_index_cohort.R` — Extract patients with index disease, assign
  index date, split baseline/follow-up
- `src/build_cohort_matrix.jl` — Build sparse feature matrix (COO) for cohort,
  applying preprocessing filters
- `src/run_onlinepca.jl` — Exact out-of-core PCA via OnlinePCA.jl
- `src/prepare_landscaper_input.jl` — Binarize PCA scores to ±1 for Landscaper
- `src/build_patient_scores.R` — Build patient_scores.tsv (Scores + Landscaper
  output + tekiyo + hospitalization, one row per patient)
- `src/cox_sanity_check_simple.R` — Cox: mortality/hospitalization ~ energy +
  age + sex
- `src/plot_pca_pairs_disease.jl` — Per-disease PCA pair density plots (Julia/
  CairoMakie, 11 PNGs per disease)
- `src/cox_per_basin.R` — Per-basin Cox: energy ~ mortality within each basin
- `src/km_per_basin.R` — KM survival curves per basin, energy tertile coloring
- `src/tagcloud_per_basin.R` — Tag cloud of enriched comorbidity codes per basin
  (reads binary COO: header nrows/ncols as Int32, then row_id/col_id pairs;
  English names from data/who_icd10_en.tsv)
- `src/tagcloud_per_state.R` — Tag cloud per Allstates pattern (up to 128),
  same COO/English name approach as tagcloud_per_basin.R
- `src/plot_pca_pairs_fullperiod.jl` — Full-period PCA pair density plots
  (reads Scores.csv directly, all observation rows)
- `src/plot_landscape_atlas.R` — Landscape visualization, KM curves, forest
  plots, atlas heatmap
- `src/trajectory_to_death.R` — Trajectory-to-death analysis: state sequence,
  energy, basin membership for deceased patients. Heatmap of state occupancy
  by distance-to-death. Requires preprocessed window lookup
  (`data/rowid_window_lookup.tsv`)
- `src/build_window_lookup.sh` — Preprocess row_id_number_small.txt into
  row_id→window TSV for trajectory analysis

## Container Images

- `ghcr.io/rikenbit/onlinepcajl:df65589` — OnlinePCA.jl
- `ghcr.io/chiba-ai-med/landscaper:main` — Landscaper
- `docker://koki/mmjl:20250725` — Julia + MatrixMarket
- `docker://koki/desc_onlinexxx_landscaper:20260602` — Julia + R analysis

## Output Paths

```
output/{disease}/{preproc}/cohort/          — index_cohort.tsv, feature matrix
output/{disease}/{preproc}/pca/{dim}/       — Scores.mm, Eigen_vectors.csv
output/{disease}/{preproc}/landscaper/{dim}/ — Allstates, Basin, Energy
output/{disease}/{preproc}/cox/{dim}/       — energy_sanity_check.tsv,
                                             energy_sanity_check_per_basin.tsv
output/{disease}/{preproc}/patient/         — patient_scores.tsv
summary/atlas_summary.tsv                  — cross-setting summary
plot/{disease}/pairs_*.png                 — PCA pair density plots
plot/{disease}/km_per_basin.png            — KM curves per basin
plot/{disease}/tagcloud_basin*.png         — tag cloud per basin (English names)
plot/{disease}/basin_code_enrichment.tsv   — enrichment table per basin
plot/{disease}/tagcloud_state*.png         — tag cloud per Allstates pattern
plot/{disease}/state_code_enrichment.tsv   — enrichment table per state
plot/{disease}/pairs_fullperiod.png        — PCA pair density (all periods)
plot/{disease}/trajectory_heatmap.png     — state occupancy by distance-to-death
plot/{disease}/trajectory_energy.png      — energy trajectory before death
plot/{disease}/trajectory_basin_frac.png  — fraction in basin before death
plot/{disease}/trajectories_dead.tsv      — per-patient trajectory data
plot/{disease}/{preproc}/{dim}/            — landscape figures
plot/atlas_heatmap.png                     — PASS/FAIL heatmap across diseases
data/rowid_window_lookup.tsv              — preprocessed row_id→window mapping
```

## Binary COO Format

`feature_matrix.coo.zst` is zstd-compressed binary:
- Header: 2 × Int32 LE (nrows, ncols)
- Body: pairs of Int32 LE (row_id, col_id), 1-indexed
- `col_ids.txt` maps `code,col_id` (remapped 1..ncols)
- `row_ids.txt` maps `old_row_id,new_row_id` (old = upstream COO row)

## Reference

The old repository `../DeSC-OnlineXXX-Landscaper` contains reusable
implementations for COO/zst I/O, OnlinePCA execution, Landscaper shell
wrappers, Cox models, and plotting functions.
