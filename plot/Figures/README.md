# Paper Figures — Disease Landscape Atlas (Patterns)

Updated: 2026-11-02

## Figure Structure (8 figures + Table 1)

| Fig/Table | Disease/Topic | Category | Key message |
|-----------|--------------|----------|-------------|
| 1 | Workflow | — | ICD-10 → PCA → binarize → 128 states → basins → Cox |
| 2 | PCA (all 51) | — | PC1=age, PC7=sex; 7 PCs capture median 12.5% of total variance |
| 3 | Atlas overview | — | Majority criterion satisfied in 27/51 (mortality), 35/51 (hospitalization) |
| 4 | Jaccard clustering | — | Basin code similarity clusters diseases into clinically valid groups |
| 5 | Trajectory | — | Basins are persistent (98% first-to-last match in 1-basin); non-basin states transient |
| Table 1 | C1–C6 classification | All | 51 diseases classified into 6 behavioral categories |
| 6 | Diabetes | C5 | 4 basins; cardiovascular vs dental; high transition asymmetry |
| 7 | Dementia + comparison | C3 | 2 basins; divergent mortality; transition comparison C3/C4/C5 |
| 8 | Tuberculosis | C2 | Single basin; reversed direction; PC6 as sole discriminant |

## Narrative flow

```
Fig 1-2  Method (workflow + PCA)
Fig 3    Atlas overview (全体俯瞰: energy-mortality scatter)
Fig 4    Jaccard clustering (外的妥当性: basin code類似度で疾患クラスタ再現)
Fig 5    Trajectory (Basin動態)
Table 1  C1–C6 classification (3軸 × 51疾患)
  → C1 (10): single basin, majority criterion satisfied — text only
  → C2 (7): single basin, majority criterion not satisfied
Fig 8      Tuberculosis (C2 representative: sankey + PC sign pattern)
  → C3 (2): multi-basin, low asymmetry
  → C4 (6): multi-basin, moderate asymmetry
  → C5 (9): multi-basin, high asymmetry
Fig 6      Diabetes (C5 representative: logitdiff + KM + discon + sankey)
Fig 7      Dementia (C3 representative: logitdiff + KM + transition comparison C3/C4/C5)
  → C6 (17): multi-basin, majority criterion not satisfied
```

---

## Figure Files

### Fig 1: Overview of the DeSC Disease Landscape Atlas analysis workflow
- `Fig01/Fig01_workflow.png`

### Fig 2: PCA captures age and sex axes across disease families
- `Fig02/Fig02a_diabetes_PC1_PC2_age.png` — (a) Diabetes PC1 vs PC2, age-colored
- `Fig02/Fig02b_uti_PC1_PC7_sex.png` — (b) UTI PC1 vs PC7, sex-colored
- `Fig02/Fig02c_cumulative_variance.png` — (c) Cumulative variance explained
- `Fig02/Fig02d_pc_correlation_vertical.png` — (d) PC-Age/Sex correlation heatmap

### Fig 3: Atlas of 51 disease families
- `Fig03/Fig03_atlas_overview.png`

### Fig 4: Disease clustering by basin code similarity
- `Fig04/Fig04_jaccard_clustering.png`

### Fig 5: Basin assignments persist among deceased patients
- `Fig05/Fig05a_trajectory_overview.png` — (a) Stable Basin / Stable Non-Basin / Moved
- `Fig05/Fig05b_trajectory_basin_start.png` — (b) Basin start → same / different / non-basin
- `Fig05/Fig05c_trajectory_nonbasin_start.png` — (c) Non-Basin start → basin / stay

### Fig 6: Diabetes energy landscape (C5)
- `Fig06/Fig06_a_diabetes_logitdiff.png` — (a) Basin 1 vs 2 log-odds bar chart
- `Fig06/Fig06_b_diabetes_km.png` — (b) Inter-basin KM
- `Fig06/Fig06_c_diabetes_discon.png` — (c) Disconnectivity graph
- `Fig06/Fig06_d_diabetes_sankey.png` — (d) Temporal basin flow

### Fig 7: Dementia with transition comparison (C3)
- `Fig07/Fig07_a_dementia_logitdiff.png` — (a) Basin comorbidity profiles
- `Fig07/Fig07_b_dementia_km.png` — (b) KM mortality + hospitalization
- `Fig07/Fig07_comparison.png` — (c) Transition comparison: Dementia(C3)/HF(C4)/Diabetes(C5)
- `Fig07/Fig07_d_dementia_sankey.png` — (d) Temporal basin flow

### Fig 8: Tuberculosis (C2)
- `Fig08/Fig08_a_tuberculosis_sankey.png` — (a) Temporal basin composition
- `Fig08/Fig08_b_pcpattern.png` — (b) PC sign pattern heatmap
- `Fig08/Fig08_c_pc6_loading.png` — (c) PC6 loading bar chart

## Supplementary Figures

See `Supplementary/README.md` for full descriptions.

- `Supplementary/FigS1_Fig03_atlas_hospitalization.png` — S1: Hospitalization atlas
- `Supplementary/FigS2_Table1_pattern_classification_space.png` — S2: Transition rate vs asymmetry
- `Supplementary/FigS3_Fig04_within_disease_basin_similarity.png` — S3: Within-disease basin similarity
- `Supplementary/FigS4_Fig04_crossdisease_discriminating_codes.png` — S4: Cross-disease codes
- `Supplementary/FigS5_sepsis_km_mortality.png` — S5a: Sepsis KM mortality
- `Supplementary/FigS5_sepsis_km_hospitalization.png` — S5b: Sepsis KM hospitalization
- `Supplementary/FigS5_sepsis_trajectory_heatmap.png` — S5c: Sepsis trajectory heatmap

## Supporting data files

- `Table1_pattern_classification.tsv` — C1–C6 classification for all 51 diseases
- `asymmetry_corrected.tsv` — Corrected asymmetry scores for multi-basin diseases
- `Supplementary/trajectory_classification.tsv` — Patient trajectory classifications
- `Supplementary/code_trajectories/` — Per-basin comorbidity code trajectory plots
