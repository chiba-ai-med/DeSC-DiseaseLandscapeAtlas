# Supplementary Figures

## FigS1. Atlas of hospitalization hazard ratios across 51 disease families

Scatter plot of observation density (x-axis) versus death rate (y-axis)
across 51 disease families, with point color indicating per-basin
majority criterion for hospitalization (HR < 1 in ≥ half of basins;
blue = satisfied, red = not satisfied) and point size reflecting
|log(HR)| from the hospitalization Cox model. Same axes as Fig. 3;
differs only in the endpoint used for point color and size
(hospitalization instead of mortality). The criterion is satisfied
in 35 of 51 families.

- File: `FigS1_Fig03_atlas_hospitalization.png`

## FigS2. Asymmetry score versus mean transition rate across ELA behavioral patterns

Scatter plot of mean per-step transition rate (x-axis) against mean
pairwise asymmetry score (y-axis) for all multi-basin diseases (basin 0
excluded, full observation period). Points are colored by the six
behavioral categories (C1–C6) defined in Table 1. Natural gaps in the
asymmetry score distribution (0.02–0.19 and 0.35–0.52) are visible,
supporting the three-level classification (C3, C4, C5).

- File: `FigS2_pattern_space.png`
- Script: `src/regen_figs2_pattern_space.py`

## FigS3. Within-disease basin-to-basin Jaccard similarity

Dot plot of pairwise Jaccard similarity between basins within each
multi-basin disease, computed from the top enriched comorbidity codes
per basin. Higher similarity indicates that basins share comorbidity
profiles; lower similarity indicates clinically distinct subgroups.
Supplements Fig. 4.

- File: `FigS3_basin_similarity.png`
- Script: `src/regen_figs3_within_disease_jaccard.R`

## FigS4. Cross-disease basin-discriminating comorbidity codes

Heatmap of the top 10 comorbidity codes with the highest mean basin
spread (difference in prevalence between the most- and least-enriched
basins) across 24 multi-basin diseases. Values are spread in percentage
points. Codes such as unspecified diabetes mellitus and chronic
periodontitis discriminate basins across many diseases, suggesting
shared comorbidity axes. Supplements Fig. 4.

- File: `FigS4_Fig04_crossdisease_discriminating_codes.png`

## FigS5. Sepsis — Kaplan–Meier survival curves for all-cause mortality (C6)

Kaplan–Meier survival curves for all-cause mortality in the sepsis
cohort, a representative of C6 (multi-basin, majority criterion not
satisfied), stratified by basin membership (Basin 1–3; n = 132,559).

- File: `FigS5_sepsis_mortality.png`

## FigS6. Sepsis — Kaplan–Meier curves for disease-specific hospitalization (C6)

Kaplan–Meier curves for disease-specific hospitalization in the sepsis
cohort, stratified by basin membership (Basin 1–3; n = 369,717,
excluding pre-landmark admissions).

- File: `FigS6_sepsis_hospitalization.png`

## FigS7. Sepsis — state occupancy heatmap before death (C6)

State occupancy heatmap before death in the sepsis cohort. Rows are
the 25 most frequent Allstates patterns among deceased patients
(n = 61,074), sorted by peak timing (states that predominate earlier
in the trajectory appear at top; states concentrated near death appear
at bottom). Star marks basin states.

- File: `FigS7_sepsis_trajectory.png`

## Supplementary Data

- `code_trajectories/` — Per-basin comorbidity code trajectory plots
  (prevalence over time for individual ICD-10 codes within each basin)
  for selected diseases (anemia, dementia, thyroid disorder, tuberculosis).
- `trajectory_classification.tsv` — Classification of individual patient
  trajectories (stable in basin, stable in non-basin, moved) for all
  multi-basin diseases.
