#!/bin/bash
# Re-run Landscaper v1.7.1 (basin fix) + downstream for all diseases
# Usage: bash src/run_landscaper_rerun.sh

set -e

CONTAINER_LS="docker://ghcr.io/chiba-ai-med/landscaper:main"
CONTAINER_R="docker://koki/desc_onlinexxx_landscaper:20260602"
TEKIYO="data/tekiyo.csv"
HOSP="data/person_hospitalization.tsv"
PREPROC="all_comorbid"
DIM=7

log() { echo "[$(date +%H:%M:%S)] $1"; }

for d in output/*/; do
    fid=$(basename "$d")
    base="output/${fid}/${PREPROC}"
    scores_mm="${base}/pca/${DIM}/Scores.mm"
    ls_dir="${base}/landscaper/${DIM}"

    if [ ! -f "$scores_mm" ]; then
        log "SKIP $fid (no Scores.mm)"
        continue
    fi

    # Phase 1: Re-run Landscaper
    log "$fid: Phase 1 - Landscaper re-run"
    tmpdir=$(mktemp -d "./landscaper_tmp.XXXXXX")
    singularity exec "$CONTAINER_LS" bash -c "
        cd '$tmpdir'
        cp /landscaper .
        cp /Snakefile .
        cp -rf /src .
        ./landscaper -i '../${scores_mm}' -o '../${ls_dir}' \
            --input_sparse=TRUE --cores=1 --memgb=100
    " 2>&1 | tail -5
    rm -rf "$tmpdir"
    log "$fid: Landscaper done"

    # Phase 2: Rebuild patient_scores.tsv
    log "$fid: Phase 2 - rebuild patient_scores.tsv"
    singularity exec "$CONTAINER_R" \
        Rscript src/build_patient_scores.R \
        "$fid" "$base" "$TEKIYO" "$HOSP" \
        "${base}/patient/patient_scores.tsv" 2>&1 | tail -5
    log "$fid: patient_scores done"

    # Phase 3: Re-run Cox energy sanity check
    log "$fid: Phase 3 - Cox sanity check"
    mkdir -p "${base}/cox/${DIM}"
    singularity exec "$CONTAINER_R" \
        Rscript src/cox_energy_sanity_check.R \
        "$fid" "$DIM" "$base" "$TEKIYO" "${base}/cox/${DIM}" 2>&1 | tail -3
    log "$fid: Cox done"

    log "$fid: ALL PHASES COMPLETE"
    echo ""
done

# Phase 4: Global summaries
log "Phase 4 - Basin PASS/FAIL profile"
singularity exec "$CONTAINER_R" \
    Rscript src/basin_pass_fail_profile.R 2>&1 | tail -5

log "Phase 5 - Compile atlas summary"
singularity exec "$CONTAINER_R" \
    Rscript src/compile_summary.R 2>&1 | tail -5

# Phase 6: Re-draw basin pair plots (all diseases, single Julia session)
log "Phase 6 - Re-draw basin pair plots"
singularity exec "$CONTAINER_R" bash src/run_all_pca_pairs_basin_only.sh 2>&1 | tail -10

log "ALL DONE"
