#!/bin/bash
set -euo pipefail

# Full pipeline: 51 disease families, Phase 1-5
# Processes smallest-to-largest for faster initial results.
# Skips already-completed steps.

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$BASEDIR"

COO_ZST="data/coo.txt.zst"
COL_IDS="data/col_id_disease_name_small.txt"
FAMILIES_TSV="config/disease_families.tsv"
ROW_IDS_ORIG="data/row_id_number_small.txt"
TEKIYO="data/tekiyo.csv"
HOSP="data/person_hospitalization.tsv"

CONTAINER_PCA="docker://ghcr.io/rikenbit/onlinepcajl:df65589"
CONTAINER_MMJL="docker://koki/mmjl:20250725"
CONTAINER_LANDSCAPER="docker://ghcr.io/chiba-ai-med/landscaper:main"
CONTAINER_R="docker://koki/desc_onlinexxx_landscaper:20260602"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

###############################################################
# Phase 1+2: Extract all cohort COOs
###############################################################
phase12_extract() {
  # Check if all COOs exist
  local all_done=true
  while IFS=$'\t' read -r fid fname prefixes rest; do
    [ "$fid" = "family_id" ] && continue
    if [ ! -f "output/${fid}/all_comorbid/cohort/feature_matrix.coo.zst" ]; then
      all_done=false
      break
    fi
  done < "$FAMILIES_TSV"

  if $all_done; then
    log "Phase 1+2: All COOs already extracted, skipping"
    return
  fi

  log "Phase 1+2: Extracting cohort COOs for all families..."
  python3 src/extract_cohort_coo.py \
    "$COO_ZST" "$COL_IDS" "$FAMILIES_TSV" output/
  log "Phase 1+2: Done"
}

###############################################################
# Phase 1.5: Build person_map for all families
###############################################################
phase15_person_map() {
  local all_done=true
  while IFS=$'\t' read -r fid fname prefixes rest; do
    [ "$fid" = "family_id" ] && continue
    if [ ! -f "output/${fid}/all_comorbid/cohort/person_map.tsv" ]; then
      all_done=false
      break
    fi
  done < "$FAMILIES_TSV"

  if $all_done; then
    log "Phase 1.5: All person_maps exist, skipping"
    return
  fi

  log "Phase 1.5: Building person_map for all families..."
  python3 src/rebuild_row_ids.py \
    "$COO_ZST" "$COL_IDS" "$FAMILIES_TSV" "$ROW_IDS_ORIG" output/
  log "Phase 1.5: Done"
}

###############################################################
# Phase 3-5 per family (sorted by size)
###############################################################
process_family() {
  local fid="$1"
  local base="output/${fid}/all_comorbid"

  log "--- $fid: starting ---"

  # Phase 3a: ooc_cov
  if [ ! -f "$base/ooc_cov/cov_mat.csv" ]; then
    log "$fid: Phase 3a - ooc_cov"
    mkdir -p "$base/ooc_cov"
    singularity exec "$CONTAINER_PCA" bash -c "
      export JULIA_DEPOT_PATH='/usr/local/julia'
      export JULIA_HISTORY='/dev/null'
      julia src/run_ooc_cov.jl \
        $base/cohort/feature_matrix.coo.zst \
        $base/ooc_cov/cov_mat.csv \
        $base/ooc_cov/colmeanvec.csv
    " 2>&1 | grep -v "warn rootless" | tail -3
  else
    log "$fid: ooc_cov exists, skipping"
  fi

  # Phase 3b: eigen
  if [ ! -f "$base/ooc_cov/eigenvalues.csv" ]; then
    log "$fid: Phase 3b - eigen"
    singularity exec "$CONTAINER_R" \
      Rscript src/eigen.R \
        "$base/ooc_cov/cov_mat.csv" \
        "$base/ooc_cov/eigenvalues.csv" \
        "$base/ooc_cov/eigenvectors.csv" 2>&1 | tail -2
  else
    log "$fid: eigen exists, skipping"
  fi

  # Phase 3c: projection (dim=7 only)
  if [ ! -f "$base/pca/7/Scores.mm" ]; then
    log "$fid: Phase 3c - projection dim=7"
    # Decompress COO
    zstd -d -f -k "$base/cohort/feature_matrix.coo.zst" \
      -o "$base/cohort/feature_matrix.coo"
    mkdir -p "$base/pca/7"
    singularity exec "$CONTAINER_MMJL" bash -c "
      mkdir -p /tmp/julia_depot
      export JULIA_DEPOT_PATH='/tmp/julia_depot:/julia_packages'
      export JULIA_HISTORY='/dev/null'
      /usr/local/julia/bin/julia src/project_scores.jl \
        $base/cohort/feature_matrix.coo \
        $base/ooc_cov/eigenvectors.csv \
        $base/ooc_cov/colmeanvec.csv \
        $base/pca \
        7
    " 2>&1 | grep -v "warn rootless" | tail -3
    # Clean up decompressed COO
    rm -f "$base/cohort/feature_matrix.coo"
  else
    log "$fid: projection exists, skipping"
  fi

  # Phase 4: Landscaper dim=7
  if [ ! -f "$base/landscaper/7/Allstates.tsv" ]; then
    log "$fid: Phase 4 - Landscaper dim=7"
    mkdir -p "$base/landscaper/7"
    local tmpdir
    tmpdir=$(mktemp -d "./landscaper_tmp.XXXXXX")
    singularity exec "$CONTAINER_LANDSCAPER" bash -c "
      cd '$tmpdir'
      cp /landscaper .
      cp /Snakefile .
      cp -rf /src .
      ./landscaper -i '../$base/pca/7/Scores.mm' -o '../$base/landscaper/7' \
        --input_sparse=TRUE --cores=1 --memgb=100
    " 2>&1 | tail -3
    rm -rf "$tmpdir"
  else
    log "$fid: Landscaper exists, skipping"
  fi

  # Phase 5a: Cox sanity check
  if [ ! -f "$base/cox/7/energy_sanity_check.tsv" ]; then
    log "$fid: Phase 5a - Cox sanity check"
    mkdir -p "$base/cox/7"
    singularity exec "$CONTAINER_R" \
      Rscript src/cox_sanity_check_simple.R \
        "$fid" 7 "$base" "$TEKIYO" "$base/cox/7" 2>&1 | tail -5
  else
    log "$fid: Cox exists, skipping"
  fi

  # Phase 5b: Build patient scores
  if [ ! -f "$base/patient/patient_scores.tsv" ]; then
    log "$fid: Phase 5b - patient scores"
    mkdir -p "$base/patient"
    singularity exec "$CONTAINER_R" \
      Rscript src/build_patient_scores.R \
        "$fid" "$base" "$TEKIYO" "$HOSP" \
        "$base/patient/patient_scores.tsv" 2>&1 | tail -3
  else
    log "$fid: patient scores exists, skipping"
  fi

  # Phase 5c: Plots
  if [ ! -f "plot/${fid}/pairs_${fid}_death.png" ]; then
    log "$fid: Phase 5c - PCA pairs plots"
    mkdir -p "plot/${fid}"
    singularity exec "$CONTAINER_R" \
      Rscript src/plot_pca_pairs.R \
        "$base/patient/patient_scores.tsv" "$fid" "plot/${fid}" 2>&1 | tail -3
  else
    log "$fid: PCA plots exist, skipping"
  fi

  if [ ! -f "plot/${fid}/resid_energy_death_${fid}.png" ]; then
    log "$fid: Phase 5d - surprise analysis"
    singularity exec "$CONTAINER_R" \
      Rscript src/plot_surprises.R \
        "$base/patient/patient_scores.tsv" \
        "$base/ooc_cov/eigenvectors.csv" \
        "$base/cohort/col_ids.txt" \
        "$fid" "plot/${fid}" 2>&1 | tail -5
  else
    log "$fid: surprise analysis exists, skipping"
  fi

  log "--- $fid: complete ---"
}

###############################################################
# Main
###############################################################
log "=========================================="
log "Full pipeline: 51 disease families"
log "=========================================="

# Phase 1+2
phase12_extract

# Phase 1.5
phase15_person_map

# Get family sizes and sort smallest-first
log "Sorting families by cohort size..."
declare -A FAMILY_SIZES
while IFS=$'\t' read -r fid fname prefixes rest; do
  [ "$fid" = "family_id" ] && continue
  stats="output/${fid}/all_comorbid/cohort/matrix_stats.tsv"
  if [ -f "$stats" ]; then
    n=$(grep "^n_patients" "$stats" | cut -f2)
    FAMILY_SIZES[$fid]=$n
  else
    FAMILY_SIZES[$fid]=0
  fi
done < "$FAMILIES_TSV"

# Sort by size
SORTED_FAMILIES=$(for fid in "${!FAMILY_SIZES[@]}"; do
  echo "${FAMILY_SIZES[$fid]} $fid"
done | sort -n | awk '{print $2}')

log "Processing order (smallest first):"
for fid in $SORTED_FAMILIES; do
  log "  ${fid}: ${FAMILY_SIZES[$fid]} rows"
done

# Phase 3-5 for each family
for fid in $SORTED_FAMILIES; do
  process_family "$fid"
done

# Final summary
log "=========================================="
log "Compiling final summary"
log "=========================================="
mkdir -p summary
singularity exec "$CONTAINER_R" \
  Rscript src/compile_summary.R summary/atlas_summary.tsv 2>&1 | tail -10

log "=========================================="
log "ALL DONE"
log "=========================================="
