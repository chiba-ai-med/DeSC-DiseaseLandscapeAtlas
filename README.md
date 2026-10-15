# DeSC-DiseaseLandscapeAtlas

Disease-family-specific Energy Landscape Analysis (ELA) on DeSC claims data.

## Motivation

The prior whole-cohort ELA (`DeSC-OnlineXXX-Landscaper`) produced a landscape
whose low-energy basins reflected high-frequency lifestyle diseases. The
resulting Energy Cox showed *higher* energy = *higher* risk — essentially a
"recording atypicality" score rather than a disease-progression landscape.

This project rebuilds the analysis with **disease-family-specific index
cohorts** and explores preprocessing strategies that yield landscapes where
**low-energy basins correspond to clinically severe / post-complication states**
(i.e., Energy HR < 1 in the mortality/hospitalization Cox).

## Pilot Disease Families

| Key        | Name               | ICD-10 Codes             |
|------------|--------------------|--------------------------|
| dementia   | Dementia           | F00-F03, G30             |
| diabetes   | Diabetes mellitus  | E10-E14                  |
| hf_cancer  | Heart failure / Cancer | I50, C00-C97         |

## Preprocessing Strategies

| Setting       | Description                                        |
|---------------|----------------------------------------------------|
| all_comorbid  | All comorbidity codes (index disease excluded)      |
| drop_top10    | Remove top-10 most frequent codes in cohort         |
| drop_top25    | Remove top-25 most frequent codes in cohort         |
| drop_top50    | Remove top-50 most frequent codes in cohort         |
| drop_top100   | Remove top-100 most frequent codes in cohort        |
| complications | Restrict to severity/complication-related codes     |

## Workflow

```bash
# 1. Define index cohorts
snakemake -s workflow/cohort.smk --cores <n>

# 2. Build feature matrices
snakemake -s workflow/matrix.smk --cores <n>

# 3. Run OnlinePCA
snakemake -s workflow/onlinepca.smk --cores <n>

# 4. Run Landscaper (ELA)
snakemake -s workflow/landscaper.smk --cores <n>

# 5. Cox energy sanity check
snakemake -s workflow/cox.smk --cores <n>

# 6. Compile summary table
snakemake -s workflow/summary.smk --cores <n>

# 7. Generate figures
snakemake -s workflow/figures.smk --cores <n>
```

## Output Structure

```
output/{disease_family}/{preproc}/
  cohort/          # index cohort + feature matrix
  pca/{dim}/       # PCA scores, loadings
  landscaper/{dim}/ # Allstates, Basin, SubGraph, Energy
  cox/{dim}/       # Energy sanity check results
summary/
  atlas_summary.tsv  # Cross-disease x preprocessing summary
plot/
  {disease_family}/{preproc}/{dim}/  # Landscape + Cox figures
```

## Requirements

- Snakemake >= 8.11.1
- Julia (OnlinePCA.jl)
- R (survival, data.table, ggplot2)
- Landscaper container: `ghcr.io/chiba-ai-med/landscaper:main`
