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
![Figure 3]('"$FIGS"'/Fig03_atlas_overview.png)\
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
![Figure 6]('"$FIGS"'/Fig06/Fig06_thyroid.png)\
' \
  -e '/^### Figure 7\./i\
![Figure 7]('"$FIGS"'/Fig07/Fig07a_dementia_paired_km.png)\
' \
  -e '/^### Figure 8\./i\
![Figure 8]('"$FIGS"'/Fig08/Fig08_anemia.png)\
' \
  -e '/^### Figure 9\./i\
![Figure 9]('"$FIGS"'/Fig09/Fig09_tuberculosis.png)\
' \
  -e '/^### Figure S1\./i\
![Figure S1]('"$FIGS"'/FigS_atlas_hospitalization.png)\
' \
  -e '/^### Figure S2\./i\
![Figure S2]('"$FIGS"'/Supplementary/FigS_jaccard_2blocks.png)\
' \
  -e '/^### Figure S3\./i\
![Figure S3]('"$FIGS"'/Supplementary/FigS_jaccard_basin_level.png)\
' \
  -e '/^### Figure S4\./i\
![Figure S4]('"$FIGS"'/Supplementary/FigS_pattern_space.png)\
' \
  -e '/^### Figure S5\./i\
![Figure S5]('"$FIGS"'/Supplementary/FigS_death_hosp.png)\
' \
  -e '/^### Figure S6\./i\
![Figure S6]('"$FIGS"'/Supplementary/FigS_crossdisease_codes.png)\
' \
  -e '/^### Figure S7\./i\
![Figure S7]('"$FIGS"'/Supplementary/Fig06_diabetes.png)\
' \
  -e '/^### Figure S8\./i\
![Figure S8a: Sepsis KM mortality]('"$FIGS"'/Supplementary/FigS8_sepsis_km_mortality.png)\
\
![Figure S8b: Sepsis KM hosp]('"$FIGS"'/Supplementary/FigS8_sepsis_km_hosp.png)\
\
![Figure S8c: Sepsis tagcloud B1]('"$FIGS"'/Supplementary/FigS8_sepsis_tagcloud_B1.png)\
\
![Figure S8d: Sepsis tagcloud B2]('"$FIGS"'/Supplementary/FigS8_sepsis_tagcloud_B2.png)\
\
![Figure S8e: Sepsis tagcloud B3]('"$FIGS"'/Supplementary/FigS8_sepsis_tagcloud_B3.png)\
\
![Figure S8f: Sepsis trajectory]('"$FIGS"'/Supplementary/FigS8_sepsis_traj_heatmap.png)\
' \
  -e '/^### Figure S9\./i\
![Figure S9a: Dementia B1 death]('"$FIGS"'/Supplementary/FigS_code_traj_dementia_B1.png)\
\
![Figure S9b: Dementia B2 death]('"$FIGS"'/Supplementary/FigS_code_traj_dementia_B2.png)\
\
![Figure S9c: Thyroid B1 death]('"$FIGS"'/Supplementary/FigS_code_traj_thyroid_disorder_B1.png)\
\
![Figure S9d: Thyroid B2 death]('"$FIGS"'/Supplementary/FigS_code_traj_thyroid_disorder_B2.png)\
\
![Figure S9e: Thyroid B3 death]('"$FIGS"'/Supplementary/FigS_code_traj_thyroid_disorder_B3.png)\
\
![Figure S9f: Thyroid B4 death]('"$FIGS"'/Supplementary/FigS_code_traj_thyroid_disorder_B4.png)\
\
![Figure S9g: Anemia B1 death]('"$FIGS"'/Supplementary/FigS_code_traj_anemia_B1.png)\
\
![Figure S9h: Anemia B3 death]('"$FIGS"'/Supplementary/FigS_code_traj_anemia_B3.png)\
\
![Figure S9i: TB B1 death]('"$FIGS"'/Supplementary/FigS_code_traj_tuberculosis_B1.png)\
\
![Figure S9j: Dementia B1 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_dementia_B1.png)\
\
![Figure S9k: Dementia B2 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_dementia_B2.png)\
\
![Figure S9l: Thyroid B1 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_thyroid_disorder_B1.png)\
\
![Figure S9m: Thyroid B2 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_thyroid_disorder_B2.png)\
\
![Figure S9n: Thyroid B3 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_thyroid_disorder_B3.png)\
\
![Figure S9o: Thyroid B4 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_thyroid_disorder_B4.png)\
\
![Figure S9p: Anemia B1 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_anemia_B1.png)\
\
![Figure S9q: Anemia B3 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_anemia_B3.png)\
\
![Figure S9r: TB B1 hosp]('"$FIGS"'/Supplementary/FigS_code_traj_hosp_tuberculosis_B1.png)\
' \
  manuscript.md > "$TMP_MD"

pandoc "$TMP_MD" -o manuscript.docx --standalone 2>&1
rm -f "$TMP_MD"
echo "Done: manuscript.docx ($(du -h manuscript.docx | cut -f1))"
