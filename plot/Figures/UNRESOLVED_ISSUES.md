# Unresolved Issues — Figure/Table Deliverables

Status as of 2026-11-02.

## Completed

1. **Pattern names renamed to C1–C6** in manuscript text, Table 1, figure legends, supplementary legends, Discussion
2. **Case study figures realigned**: Fig 6=Diabetes (C5), Fig 7=Dementia (C3), Fig 8=Tuberculosis (C2), Fig 9 removed
3. **Supplementary legends rewritten (S1–S5)**: FigS1 axes/color, FigS5 basin list and n values corrected
4. **build_docx.sh updated**: composite images for Fig 2 (Fig02_pca.png) and Fig 5 (Fig05_trajectory.png)
5. **Text inconsistencies fixed**: "net directional flow" removed, full-period vs last-20 clarified, Methods §4.7 trajectory storage distinguished from visualization
6. **Diabetes numbers synchronized**: §2.6 body and Fig 6 legend now use per_basin_mortality_summary.tsv (landmark design) HR and death rates, with deaths/N shown explicitly. Fig 6a n values identified as full-cohort (patient_scores.tsv) and distinguished in legend.
7. **FigS3 Jaccard values verified**: All 56 within-disease pairs independently recomputed from basin_code_enrichment.tsv (top-30 codes by freq_in) and matched jaccard_all_basins.tsv exactly (0 mismatches). Pair counts match C(n_distinct(basin_id>0), 2) for every disease.

## Author action items — specific changes needed

### 8. Fig 6b: replace KM with landmark design

Current `Fig06/Fig06_b_diabetes_km.png` uses patient_scores.tsv with total insurance period (observable_end_ym − observable_start_ym) as survival time. Manuscript describes landmark design (fu_time from landmark window end).

- **Candidate**: `paper/claude/review_candidates/Fig06_b_diabetes_km_landmark.png`
- **Generation code**: `paper/claude/review_candidates/fig06b_landmark_km.R`
- **Input**: `output/diabetes/all_comorbid/landmark/landmark_cohort.tsv` (existing file, no recalculation)
- **Population**: shibou_flg_riyouka=1, basin_id>0, fu_time>0 → n=1,839,848
- **Decision needed**: adopt the landmark candidate or keep the current figure

### 9. Fig 6a: Basin 2 HR label correction

Image text shows Basin 2 HR = 1.17. Per-basin landmark Cox result (per_basin_mortality_summary.tsv) is HR = 1.193.

- **File**: `Fig06/Fig06_a_diabetes_logitdiff.png`
- **Location**: Basin 2 panel, HR annotation
- **Change**: 1.17 → 1.19
- Basin 1 shows 0.84, consistent with landmark HR = 0.844 — no change needed

### 10. FigS3: basin count labels

Y-axis labels show max(basin_id) instead of the number of distinct basins (basin_id > 0). Jaccard dot counts and values are correct (verified); only the labels are wrong.

**Root cause**: Generation code is not in the repository (image was added in commit 8c4c0a1 without a script). The observed labeling pattern — every affected disease has non-contiguous basin IDs where max > n_distinct — is consistent with max(basin_id) being used for the label, but this is an inference from symptoms, not a confirmed code reading.

**Affected diseases** (confirmed from data):

| Disease | Image label | Actual distinct basins | Basin IDs |
|---|---|---|---|
| CKD | 6 basins | 2 | {1, 6} |
| UTI | 6 basins | 3 | {1, 4, 6} |
| BPH | 4 basins | 3 | {1, 2, 4} |
| Depression | 3 basins | 2 | {1, 3} |
| Prostate cancer | 3 basins | 2 | {1, 3} |
| Anemia | 3 basins | 2 | {1, 3} |

Diseases with contiguous basin IDs (e.g., diabetes {1,2,3,4}) display correct counts.

**Jaccard verification**: 56 within-disease pairs recomputed from `plot/{disease}/basin_code_enrichment.tsv` (top-30 codes by freq_in descending) using Jaccard = |intersection| / |union|. All values matched `summary/jaccard_all_basins.tsv` exactly (tolerance < 1e-10). Pair counts per disease equal C(n_distinct, 2). **Values are correct; only the y-axis labels need correction.**

- **File**: `Supplementary/FigS3_Fig04_within_disease_basin_similarity.png`
- **Changes**: CKD 6→2, UTI 6→3, BPH 4→3, Depression 3→2, Prostate cancer 3→2, Anemia 3→2

### 11. FigS2: old pattern names → C1–C6

- **File**: `Supplementary/FigS2_Table1_pattern_classification_space.png`
- **Locations**: legend (4 entries) and right-margin region labels
- **Changes**:
  - Symmetric → C3
  - Convergent → C4
  - Progressive → C5
  - Structured reversed → C6

### 12. Fig07_comparison.png: old pattern names → C1–C6

- **File**: `Fig07/Fig07_comparison.png`
- **Locations**: panel titles (3 panels)
- **Changes**:
  - "Dementia (Symmetric)" → "Dementia (C3)"
  - "Heart Failure (Convergent)" → "Heart Failure (C4)"
  - "Diabetes (Progressive)" → "Diabetes (C5)"

No other figure images contain old pattern names (confirmed by visual inspection of all 65 PNGs).

### 13. Fig01 workflow numbers

Numbers 27/51, 35/51, 1–4 basins confirmed from data — no change needed. Author should verify these appear correctly on the workflow diagram image.

## Notes

- Fig 6/7/8 composites: author composes in Keynote from individual panels in Fig06/, Fig07/, Fig08/
- Current manuscript.docx contains updated text but embeds the pre-correction figure images. It is a review draft, not a submission-ready build. After author completes items 8–12, the Word should be rebuilt.
