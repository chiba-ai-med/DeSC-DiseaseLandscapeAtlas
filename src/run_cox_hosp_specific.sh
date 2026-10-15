#!/bin/bash
set -euo pipefail

CONTAINER="docker://koki/desc_onlinexxx_landscaper:20260602"
SCRIPT="src/cox_hospitalization_specific.R"
TEKIYO="data/tekiyo.csv"

for fid in $(ls -1 output/); do
  scores="output/${fid}/all_comorbid/patient/patient_scores.tsv"
  hosp="output/${fid}/all_comorbid/cohort/person_hospitalization_specific.tsv"
  outdir="output/${fid}/all_comorbid/cox/7"
  out_file="${outdir}/energy_sanity_check_hosp_specific.tsv"

  if [ ! -f "$scores" ]; then
    echo "SKIP $fid: no patient_scores.tsv"
    continue
  fi

  if [ ! -f "$hosp" ]; then
    echo "SKIP $fid: no person_hospitalization_specific.tsv"
    continue
  fi

  if [ -f "$out_file" ]; then
    echo "SKIP $fid: already done"
    continue
  fi

  echo "=== $fid ==="
  mkdir -p "$outdir"
  singularity exec "$CONTAINER" Rscript "$SCRIPT" "$fid" "$scores" "$hosp" "$TEKIYO" "$out_file"
done

echo "=== ALL DONE ==="
