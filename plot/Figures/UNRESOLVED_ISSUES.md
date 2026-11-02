# Unresolved Issues — Figure/Table Deliverables

Status as of 2026-11-02.

## Completed in this session (2026-11-02)

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

4. **build_docx.sh updated** with correct figure paths

## Open issues

### Author action needed

5. **Figure text still shows old pattern names**
   - Images cannot be modified per author rule
   - The following figures contain old pattern names in their image text:
     - FigS2: legend may show "Symmetric", "Convergent", "Progressive" etc.
     - FigS5: title may say "Structured reversed"
   - **Action**: Author should update figure images where pattern names appear

6. **FigS3 basin counts in image**
   - Image labels show CKD "6 basins" and UTI "6 basins", but actual data has 2 and 3 respectively
   - This is from old data before basin 0 exclusion fix; image cannot be changed
   - **Action**: Author should regenerate FigS3 or note in manuscript

7. **Fig01 workflow numbers**
   - Must show: 27/51 mortality, 35/51 hospitalization, 1–4 basins
   - **Action**: Author should verify numbers on the workflow diagram

### Low priority

8. **Composite figures**
   - Author composes multi-panel figures in Keynote
   - Individual panels are in Fig06/, Fig07/, Fig08/ directories
