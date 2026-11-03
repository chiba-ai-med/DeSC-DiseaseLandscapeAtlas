#!/bin/bash
# Build manuscript.docx with embedded figures
set -e
cd /home/koki/dev/DeSC-DiseaseLandscapeAtlas/paper/claude

FIGS=../../plot/Figures
TMP_MD=$(mktemp /tmp/manuscript_with_figs_XXXXXX.md)

# Insert figure images before each figure legend
# Fig 5–8: author-composed single images (not individual panels)
sed \
  -e '/^### Figure 1\./i\
![Figure 1]('"$FIGS"'/Fig01/Fig01_workflow.png)\
' \
  -e '/^### Figure 2\./i\
![Figure 2]('"$FIGS"'/Fig02/Fig02_pca.png)\
' \
  -e '/^### Figure 3\./i\
![Figure 3]('"$FIGS"'/Fig03/Fig03_atlas_overview.png)\
' \
  -e '/^### Figure 4\./i\
![Figure 4]('"$FIGS"'/Fig04/Fig04_jaccard_clustering.png)\
' \
  -e '/^### Figure 5\./i\
![Figure 5]('"$FIGS"'/Fig05/Fig05_trajectory.png)\
' \
  -e '/^### Figure 6\./i\
![Figure 6]('"$FIGS"'/Fig06/Fig06_diabetes.png)\
' \
  -e '/^### Figure 7\./i\
![Figure 7]('"$FIGS"'/Fig07/Fig07_dementia.png)\
' \
  -e '/^### Figure 8\./i\
![Figure 8]('"$FIGS"'/Fig08/Fig08.png)\
' \
  -e '/^### Figure S1\./i\
![Figure S1]('"$FIGS"'/Supplementary/FigS1_Fig03_atlas_hospitalization.png)\
' \
  -e '/^### Figure S2\./i\
![Figure S2]('"$FIGS"'/Supplementary/FigS2_pattern_space.png)\
' \
  -e '/^### Figure S3\./i\
![Figure S3]('"$FIGS"'/Supplementary/FigS3_basin_similarity.png)\
' \
  -e '/^### Figure S4\./i\
![Figure S4]('"$FIGS"'/Supplementary/FigS4_Fig04_crossdisease_discriminating_codes.png)\
' \
  -e '/^### Figure S5\./i\
![Figure S5]('"$FIGS"'/Supplementary/FigS5_sepsis_mortality.png)\
' \
  -e '/^### Figure S6\./i\
![Figure S6]('"$FIGS"'/Supplementary/FigS6_sepsis_hospitalization.png)\
' \
  -e '/^### Figure S7\./i\
![Figure S7]('"$FIGS"'/Supplementary/FigS7_sepsis_trajectory.png)\
' \
  manuscript.md > "$TMP_MD"

pandoc "$TMP_MD" -o manuscript.docx --standalone 2>&1
rm -f "$TMP_MD"
echo "Done: manuscript.docx ($(du -h manuscript.docx | cut -f1))"
