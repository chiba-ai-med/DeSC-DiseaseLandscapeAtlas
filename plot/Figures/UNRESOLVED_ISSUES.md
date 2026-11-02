# Unresolved Issues — Figure/Table Deliverables

Status as of 2026-11-02.

## Completed

1. **Pattern names renamed to C1–C6**
   - Table 1 in manuscript: Homogeneous→C1, Reversed→C2, Symmetric→C3, Convergent→C4, Progressive→C5, Structured reversed→C6
   - Table1_pattern_classification.tsv updated
   - All body text, figure legends, supplementary legends updated
   - Discussion updated (C3/C5 labels, C6 category, removed Figure S5 reference)

2. **Case study figures realigned**
   - Fig 6: Thyroid → Diabetes (C5)
   - Fig 7: Dementia paired KM → Dementia (C3) + transition comparison (C3/C4/C5)
   - Fig 8: Anemia transition trio → Tuberculosis (C2)
   - Fig 9: Removed (TB content moved to Fig 8)

3. **Supplementary legends rewritten (S1–S5)**
   - Old S2–S9 removed; current S1–S5 match Supplementary/README.md
   - S1 described as scatter plot (not heatmap)
   - S2 axes corrected (x=transition rate, y=asymmetry)
   - S5 KM coloring described as by basin (not energy tertile)

4. **build_docx.sh updated**
   - Correct figure paths, composite images for Fig 2 (Fig02_pca.png) and Fig 5 (Fig05_trajectory.png)

5. **Text inconsistencies fixed (manuscript.md)**
   - §2.5 C4: "net directional flow" → "asymmetric conditional transition probabilities"
   - §2.6 diabetes: Added absolute count comparison (B2→B1: 22,923 vs B1→B2: 21,733)
   - §2.7: Clarified full-period vs last-20 differ only in time window, both use deceased patients
   - Fig 7 legend: Same fix
   - Methods §4.7: Distinguished full trajectory storage (up to 83 positions) from 20-position visualization
   - FigS1 legend: y-axis is death rate (not hospitalization rate)
   - FigS5 legend: Basin 1–3 only (no non-basin), added n values

6. **Fig 6b landmark KM candidate generated**
   - `paper/claude/review_candidates/Fig06_b_diabetes_km_landmark.png`
   - Uses landmark_cohort.tsv with fu_time/event_death (not patient_scores.tsv)
   - Existing Fig06_b_diabetes_km.png is non-landmark design (observable_end_ym − observable_start_ym)

7. **FigS3 basin count root cause identified** (see below)

## Open issues

### Author action needed

8. **Fig 6b design mismatch**
   - Current `Fig06_b_diabetes_km.png` uses patient_scores.tsv with total insurance period as time
   - Manuscript describes landmark design (fu_time from landmark window end)
   - Candidate replacement: `paper/claude/review_candidates/Fig06_b_diabetes_km_landmark.png`
   - **Action**: Author decides whether to adopt the landmark candidate

9. **FigS3 basin count labels (image only)**
   - **Root cause**: Generation code used `max(basin_id)` instead of `n_distinct(basin_id[basin_id > 0])`. Basin IDs are non-contiguous (e.g., CKD has basins {1,6}, max=6, distinct=2).
   - **Jaccard values are UNAFFECTED**: dot counts per disease match C(n_distinct, 2), confirming computation used correct basin IDs. Only the y-axis labels are wrong.
   - No generation script in repo (commit 8c4c0a1 added image without code).
   - Affected diseases (confirmed from data):
     - CKD: shows "6 basins", actual 2 (IDs: 1,6)
     - UTI: shows "6 basins", actual 3 (IDs: 1,4,6)
     - BPH: shows "4 basins", actual 3 (IDs: 1,2,4)
     - Depression: shows "3 basins", actual 2 (IDs: 1,3)
     - Prostate cancer: shows "3 basins", actual 2 (IDs: 1,3)
     - Anemia: shows "3 basins", actual 2 (IDs: 1,3)
     - Diseases with contiguous IDs (e.g., diabetes {1,2,3,4}) are correct.
   - **Action**: Author should regenerate FigS3 with `n_distinct` for labels

10. **Figure images with old pattern names (confirmed)**
    - `FigS2_Table1_pattern_classification_space.png`:
      - Legend: "Symmetric", "Convergent", "Progressive", "Structured reversed" → should be C3, C4, C5, C6
      - Right margin labels: same old names
    - `Fig07_comparison.png`:
      - Panel titles: "Dementia (Symmetric)", "Heart Failure (Convergent)", "Diabetes (Progressive)" → should be "Dementia (C3)", "Heart Failure (C4)", "Diabetes (C5)"
    - No other figure images contain old pattern names.
    - **Action**: Author should update these two images

11. **Fig01 workflow numbers** — confirmed 27/51, 35/51, 1–4 basins (no re-check needed)

### Low priority

12. **Composite figures for Fig 6/7/8**
    - Author composes multi-panel figures in Keynote
    - Individual panels are in Fig06/, Fig07/, Fig08/ directories
    - build_docx.sh embeds individual panels for these (no composite exists)
