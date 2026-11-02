#!/bin/bash
# Build manuscript.docx with embedded figures
set -e
cd /home/koki/dev/DeSC-DiseaseLandscapeAtlas/paper/claude

FIGS=../../plot/Figures
TMP_MD=$(mktemp /tmp/manuscript_with_figs_XXXXXX.md)

# Insert figure images before each figure legend
sed \
  -e '/^### Figure 1\./i\
![Figure 1]('"$FIGS"'/Fig01/Fig01_workflow.png)\
' \
  -e '/^### Figure 2\./i\
![Figure 2a]('"$FIGS"'/Fig02/Fig02a_diabetes_PC1_PC2_age.png)\
\
![Figure 2b]('"$FIGS"'/Fig02/Fig02b_uti_PC1_PC7_sex.png)\
\
![Figure 2c]('"$FIGS"'/Fig02/Fig02c_cumulative_variance.png)\
\
![Figure 2d]('"$FIGS"'/Fig02/Fig02d_pc_correlation_vertical.png)\
' \
  -e '/^### Figure 3\./i\
![Figure 3]('"$FIGS"'/Fig03/Fig03_atlas_overview.png)\
' \
  -e '/^### Figure 4\./i\
![Figure 4]('"$FIGS"'/Fig04/Fig04_jaccard_clustering.png)\
' \
  -e '/^### Figure 5\./i\
![Figure 5a]('"$FIGS"'/Fig05/Fig05a_trajectory_overview.png)\
\
![Figure 5b]('"$FIGS"'/Fig05/Fig05b_trajectory_basin_start.png)\
\
![Figure 5c]('"$FIGS"'/Fig05/Fig05c_trajectory_nonbasin_start.png)\
' \
  -e '/^### Figure 6\./i\
![Figure 6a]('"$FIGS"'/Fig06/Fig06_a_diabetes_logitdiff.png)\
\
![Figure 6b]('"$FIGS"'/Fig06/Fig06_b_diabetes_km.png)\
\
![Figure 6c]('"$FIGS"'/Fig06/Fig06_c_diabetes_discon.png)\
\
![Figure 6d]('"$FIGS"'/Fig06/Fig06_d_diabetes_sankey.png)\
' \
  -e '/^### Figure 7\./i\
![Figure 7a]('"$FIGS"'/Fig07/Fig07_a_dementia_logitdiff.png)\
\
![Figure 7b]('"$FIGS"'/Fig07/Fig07_b_dementia_km.png)\
\
![Figure 7c]('"$FIGS"'/Fig07/Fig07_comparison.png)\
\
![Figure 7d]('"$FIGS"'/Fig07/Fig07_d_dementia_sankey.png)\
' \
  -e '/^### Figure 8\./i\
![Figure 8a]('"$FIGS"'/Fig08/Fig08_a_tuberculosis_sankey.png)\
\
![Figure 8b]('"$FIGS"'/Fig08/Fig08_b_pcpattern.png)\
\
![Figure 8c]('"$FIGS"'/Fig08/Fig08_c_pc6_loading.png)\
' \
  -e '/^### Figure S1\./i\
![Figure S1]('"$FIGS"'/Supplementary/FigS1_Fig03_atlas_hospitalization.png)\
' \
  -e '/^### Figure S2\./i\
![Figure S2]('"$FIGS"'/Supplementary/FigS2_Table1_pattern_classification_space.png)\
' \
  -e '/^### Figure S3\./i\
![Figure S3]('"$FIGS"'/Supplementary/FigS3_Fig04_within_disease_basin_similarity.png)\
' \
  -e '/^### Figure S4\./i\
![Figure S4]('"$FIGS"'/Supplementary/FigS4_Fig04_crossdisease_discriminating_codes.png)\
' \
  -e '/^### Figure S5\./i\
![Figure S5a]('"$FIGS"'/Supplementary/FigS5_sepsis_km_mortality.png)\
\
![Figure S5b]('"$FIGS"'/Supplementary/FigS5_sepsis_km_hospitalization.png)\
\
![Figure S5c]('"$FIGS"'/Supplementary/FigS5_sepsis_trajectory_heatmap.png)\
' \
  manuscript.md > "$TMP_MD"

pandoc "$TMP_MD" -o manuscript.docx --standalone 2>&1
rm -f "$TMP_MD"
echo "Done: manuscript.docx ($(du -h manuscript.docx | cut -f1))"
