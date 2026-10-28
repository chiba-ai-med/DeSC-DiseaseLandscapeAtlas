# Unresolved Issues — Figure/Table Deliverables

Status as of 2026-10-28.

## Completed in this session

1. **FigS4 asymmetry — basin 0 contamination fixed**
   - Old: `fix_figures_v2.R` read `basin_transition_matrix.tsv` which included No basin / basin 0
   - New: `fix_figs4_asymmetry.py` reads `stepwise_transition.tsv`, excludes basin 0
   - Impact: thyroid basins 5→4, anemia 3→2, ckd 3→2; asymmetry values changed dramatically
     (e.g., ckd 0.607→0.011, dementia 0.449→0.022)
   - Output: `summary/transition_asymmetry_summary.tsv`, `plot/Figures/asymmetry_corrected.tsv`

2. **Table 1 reclassification with data-driven boundaries**
   - Old boundaries: <0.3 / 0.3-0.6 / >0.6
   - New boundaries: <0.1 / 0.1-0.45 / ≥0.45 (gap analysis of corrected values)
   - 4 diseases changed: dementia Convergent→Symmetric, hypertension Symmetric→Convergent,
     malnutrition Symmetric→Convergent, asthma Progressive→Convergent
   - Output: `plot/Figures/Table1_pattern_classification.tsv`

3. **FigS_pattern_space.png redrawn** with corrected values and new boundary lines

4. **PCA axis interpretation (50 diseases)**
   - Direction-specific and |loading|-based cross-disease analysis
   - Key findings: PC1=age/frailty (41/50), PC2=outpatient vs inpatient, PC5=metabolic/locomotor
   - Sex axis NOT fixed to any PC (distributed PC2-7, median |cor_sex|=0.180)
   - Output: 6 TSVs in `summary/pc_*.tsv`

5. **Panels copied to panels/ directories** for Fig02, Fig05, Fig06, Fig07, Fig08, Fig09

6. **Scripts saved to src/**
   - `src/fix_figs4_asymmetry.py`
   - `src/reclassify_table1.py`
   - `src/pca_crossdisease_interpretation.py`

## Open issues

### High priority (blocks submission)

7. **Fig01: numbers need verification on the workflow diagram**
   - Current image (1.2MB) was modified by `fix_fig01_definitive.py` (PIL text overlay)
   - `.bak` file (1.8MB) is the original
   - Numbers that must appear: 27/51 mortality, 35/51 hospitalization, 1-4 basins
   - **Action**: Author should verify Fig01 numbers match these values; if not,
     restore `.bak` and re-overlay correct numbers

8. **Manuscript text not yet updated** to reflect:
   - New asymmetry boundaries (Symmetric <0.1, Convergent 0.1-0.45, Progressive ≥0.45)
   - 4 diseases reclassified (dementia, hypertension, malnutrition, asthma)
   - PCA interpretation findings (PC1 age axis in 41/50 diseases, sex axis distributed)
   - **Action**: Update `paper/claude/manuscript.md` sections on Table 1,
     pattern classification, and PCA interpretation

### Medium priority (improves paper quality)

9. **PCA supplementary figure ON HOLD**
   - User paused ("ちょっとストップ") to discuss with ChatGPT
   - The visualization approach and which current Fig (Fig02?) to replace is undecided
   - Summary TSVs are ready in `summary/pc_*.tsv`

10. **Fig06a description mismatch**
    - Prompt says "Fig6aは現在、結核コードのtagcloud" but actual file is a thyroid
      basin 3 vs 4 logitdiff bar chart (from `thyroid_fig4a_final.R`)
    - No action needed — the current image is correct per the figure layout

11. **FigS_code_traj duplicate thyroid files**
    - Both `FigS_code_traj_thyroid_B*.png` (4 files) and
      `FigS_code_traj_thyroid_disorder_B*.png` (4 files) exist
    - The `thyroid` versions appear to be an older naming; `thyroid_disorder` is correct
    - **Action**: Remove `FigS_code_traj_thyroid_B*.png` duplicates

12. **Tagcloud bar charts — RESOLVED**
    - All 51 diseases now have logitdiff bar chart PNGs: `plot/{disease}/tagcloud_scored/logitdiff/barchart_logitdiff_basin{N}.png` (190 total)
    - All 51 diseases have scoring TSVs: `plot/{disease}/tagcloud_scored/logitdiff/basin_code_logitdiff.tsv`
    - Old word cloud tagclouds still exist but are NOT used in the paper

### Low priority

13. **Fig07 has 6 KM panels** — author should decide which subset to include in the
    main figure vs supplementary (paired KM, mortality landmark, hosp per-basin,
    mortality combined, hosp combined, per-basin landmark)

14. **FigS8 sepsis and Fig06 diabetes in Supplementary** — confirmed supplementary status;
    no further action unless paper structure changes

15. **Word/docx rebuild** — after manuscript text updates, run the docx build pipeline
