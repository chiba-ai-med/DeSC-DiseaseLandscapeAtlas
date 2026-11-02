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
8. **FigS3 basin count labels corrected**: Regenerated with n_distinct(basin_id > 0) for all 34 diseases. Reference line = between-disease mean (0.570). 56 dots, all pair counts verified. Script: `src/regen_figs3_within_disease_jaccard.R`.
9. **FigS2 pattern names → C3–C6**: Legend and right-margin region labels updated. 33 diseases (34 multi-basin minus schizophrenia). Script: `src/regen_figs2_pattern_space.py`.
10. **Fig07_comparison.png panel titles → C3–C5**: Dementia (C3), Heart Failure (C4), Diabetes (C5). All cell values unchanged. Script: `src/fig07_transition_comparison.R`.
11. **Fig06_a Basin 2 HR label corrected**: 1.17 → 1.19 (per_basin_mortality_summary.tsv landmark HR = 1.193). Script: `src/fig06a_diabetes_basin1_vs_2.R`.

## Author action items — specific changes needed

### 12. Fig 6b: replace KM with landmark design

Current `Fig06/Fig06_b_diabetes_km.png` uses patient_scores.tsv with total insurance period (observable_end_ym − observable_start_ym) as survival time. Manuscript describes landmark design (fu_time from landmark window end).

- **Candidate**: `paper/claude/review_candidates/Fig06_b_diabetes_km_landmark.png`
- **Generation code**: `paper/claude/review_candidates/fig06b_landmark_km.R`
- **Input**: `output/diabetes/all_comorbid/landmark/landmark_cohort.tsv` (existing file, no recalculation)
- **Population**: shibou_flg_riyouka=1, basin_id>0, fu_time>0 → n=1,839,848
- **Decision needed**: adopt the landmark candidate or keep the current figure

### 13. Fig01 workflow numbers

Numbers 27/51, 35/51, 1–4 basins confirmed from data — no change needed. Author should verify these appear correctly on the workflow diagram image.

## Notes

- Fig 6/7/8 composites: author composes multi-panel figures in Keynote from individual panels in Fig06/, Fig07/, Fig08/
- Current manuscript.docx contains updated text but embeds the pre-correction figure images. It is a review draft, not a submission-ready build. After author completes remaining items, the Word should be rebuilt.
