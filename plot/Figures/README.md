# Paper Figures — Disease Landscape Atlas (Patterns)

Updated: 2026-10-18

## Figure Structure (10 Figures)

### Fig 1: Analysis Workflow (TO BE DRAWN)
- `Fig01_workflow_PLACEHOLDER.txt`
- DeSC → 51 disease cohorts → feature matrix → PCA(7d) → binarize → ELA → Basin/Energy → downstream
- Manual creation needed (Illustrator/BioRender)

### Fig 2: PCA — Dimensionality Reduction
- `Fig02a_pca_pairs_diabetes.png` — PC1-7 pairs density plot (diabetes)
- `Fig02b_pca_pairs_thyroid.png` — PC1-7 pairs density plot (thyroid)
- TODO: cumulative variance explained across 51 diseases

### Fig 3: Energy Landscape Visualization
- `Fig03a_landscape_diabetes.png` — Landscape summary (diabetes)
- `Fig03b_km_per_basin_diabetes.png` — KM curves per basin (diabetes)
- `Fig03c_landscape_thyroid.png` — Landscape summary (thyroid)
- `Fig03d_km_per_basin_thyroid.png` — KM curves per basin (thyroid)

### Fig 4: 51-Disease Atlas Overview
- `Fig04a_atlas_heatmap.png` — PASS/FAIL heatmap across 51 diseases
- `Fig04b_atlas_3var.png` — Multi-variable atlas summary
- TODO: Cox sanity check summary panel, Basin count distribution

### Fig 5: Validation — Known Subtype Rediscovery
- `Fig05a_anemia_basin1_tagcloud.png` — Anemia Basin1: CKD/PTH codes → renal anemia
- `Fig05b_anemia_km_per_basin.png` — Anemia KM curves
- `Fig05c_hf_basin1_tagcloud.png` — Heart failure Basin1: cross-ref with CKD
- `Fig05d_ckd_basin1_tagcloud.png` — CKD Basin1: cross-ref with HF
- TODO: Jaccard similarity panel (HF×CKD = 0.33)

### Fig 6: Diabetes — 4-Basin 2-Axis Structure
- `Fig06a_diabetes_tagcloud_basin1.png` — Basin1 (metabolic: E780)
- `Fig06b_diabetes_tagcloud_basin2.png` — Basin2 (frailty: dental/eye)
- `Fig06c_diabetes_tagcloud_basin3.png` — Basin3 (metabolic: E785)
- `Fig06d_diabetes_tagcloud_basin4.png` — Basin4 (frailty: dental/eye, highest death)
- TODO: demographics summary table/bar, Ahlqvist comparison concept panel

### Fig 7: Thyroid — 5 Basins, Largest Mortality Gap
- `Fig07a_thyroid_tagcloud_basin1.png` — Basin1
- `Fig07b_thyroid_tagcloud_basin2.png` — Basin2 (highest death, frailty)
- `Fig07c_thyroid_tagcloud_basin3.png` — Basin3
- `Fig07d_thyroid_tagcloud_basin4.png` — Basin4 (lowest death, young)
- `Fig07e_thyroid_km_per_basin.png` — KM curves showing 20.4% gap

### Fig 8: Sepsis — Acute Disease Basin Structure
- `Fig08a_sepsis_tagcloud_basin1.png` — Basin1 (catheter-related infection)
- `Fig08b_sepsis_tagcloud_basin2.png` — Basin2 (typical)
- `Fig08c_sepsis_tagcloud_basin3.png` — Basin3 (ultra-elderly, 93.9% hosp)
- `Fig08d_sepsis_km_per_basin.png` — KM curves
- No BasinNA (all patients patterned)

### Fig 9: BasinNA Bidirectionality ★★★
- `Fig09_basinna_bidirectionality.png` — Forest plot: 47 diseases, HR of largest basin vs BasinNA
- Type A (red, 31 diseases): BasinNA = high risk
- Type B (blue, 16 diseases): BasinNA = low risk
- TODO: concept diagram panel, representative Type A/B comparison table

### Fig 10: Cross-Disease Frailty/Dental Signature ★★
- `Fig10a_frailty_diabetes_basin4.png` — Diabetes Basin4 tagcloud
- `Fig10b_frailty_thyroid_basin2.png` — Thyroid Basin2 tagcloud
- `Fig10c_frailty_malnutrition_basin1.png` — Malnutrition Basin1 tagcloud
- TODO: enrichment heatmap (codes × diseases), Basin2 vs Basin4 mortality comparison

## Supplementary Figures (to be organized)

Per-disease plots for all 51 diseases available in `plot/{disease}/`:
- `km_per_basin.png` — KM curves per basin
- `pairs_density.png` — PCA pairs density
- `tagcloud_basin*.png` — Tag clouds per basin
- `trajectory_*.png` — Trajectory-to-death analyses
- `all_comorbid/7/landscape_summary.png` — Landscape visualization

## Status

- [x] Fig 1: Placeholder created (needs manual drawing)
- [x] Fig 2: PCA pairs plots placed
- [x] Fig 3: Landscape + KM placed
- [x] Fig 4: Atlas heatmap placed
- [x] Fig 5: Validation tagclouds + KM placed
- [x] Fig 6: Diabetes tagclouds placed
- [x] Fig 7: Thyroid tagclouds + KM placed
- [x] Fig 8: Sepsis tagclouds + KM placed
- [x] Fig 9: BasinNA forest plot placed
- [x] Fig 10: Frailty tagclouds placed
- [ ] Multiple TODO panels need new creation (see each figure section)
