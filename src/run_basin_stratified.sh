#!/bin/bash
set -euo pipefail

CONTAINER="docker://koki/desc_onlinexxx_landscaper:20260602"
SCRIPT="src/cox_basin_stratified.R"

for fid in $(ls -1 output/); do
  scores="output/${fid}/all_comorbid/patient/patient_scores.tsv"
  outdir="plot/${fid}"
  out_check="${outdir}/cox_basin_stratified.tsv"

  if [ ! -f "$scores" ]; then
    echo "SKIP $fid: no patient_scores.tsv"
    continue
  fi

  if [ -f "$out_check" ]; then
    echo "SKIP $fid: already done"
    continue
  fi

  echo "=== $fid ==="
  singularity exec "$CONTAINER" Rscript "$SCRIPT" "$fid" "$scores" "$outdir"
done

echo "=== ALL DONE ==="
