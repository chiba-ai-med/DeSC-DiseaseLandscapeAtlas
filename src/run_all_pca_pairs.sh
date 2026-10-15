#!/bin/bash
# Run PC1-7 pair plots for all diseases
# Usage: singularity exec <container> bash src/run_all_pca_pairs.sh

export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
export JULIA_HISTORY=/dev/null
export JULIA_SCRATCH_TRACK_ACCESS=false
export JULIA_SCRATCH_SPACE=/tmp/julia_scratch
export TMPDIR=/tmp
mkdir -p /tmp/julia_depot /tmp/julia_scratch

for d in output/*/; do
    disease=$(basename "$d")
    tsv="output/${disease}/all_comorbid/patient/patient_scores.tsv"
    outdir="plot/${disease}"
    if [ ! -f "$tsv" ]; then
        echo "SKIP $disease (no patient_scores.tsv)"
        continue
    fi
    if [ -f "${outdir}/pairs_density.png" ]; then
        echo "SKIP $disease (already done)"
        continue
    fi
    echo "=== $disease ==="
    /julia_bin/bin/julia src/plot_pca_pairs_disease.jl "$disease" "$outdir"
done

echo "ALL DONE"
