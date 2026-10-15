#!/bin/bash
set -euo pipefail

# Run Landscaper + Cox sanity check for all 3 families, dim=2-7
# Usage: bash src/run_all_landscaper_cox.sh

FAMILIES="dementia diabetes hf_cancer"
DIMS="2 3 4 5 6 7"

for fam in $FAMILIES; do
  for dim in $DIMS; do
    echo "====== $fam dim=$dim ======"
    indir="output/${fam}/all_comorbid/pca/${dim}"
    outdir="output/${fam}/all_comorbid/landscaper/${dim}"
    mkdir -p "$outdir"

    if [ -f "$outdir/Allstates.tsv" ]; then
      echo "  Landscaper already done, skipping"
    else
      echo "  Running Landscaper ..."
      tmpdir=$(mktemp -d "./landscaper_tmp.XXXXXX")
      singularity exec docker://ghcr.io/chiba-ai-med/landscaper:main bash -c "
        cd '$tmpdir'
        cp /landscaper .
        cp /Snakefile .
        cp -rf /src .
        ./landscaper -i '../${indir}/Scores.mm' -o '../${outdir}' --input_sparse=TRUE --cores=1 --memgb=100
      " 2>&1 | tail -3
      rm -rf "$tmpdir"
      echo "  Landscaper done"
    fi

    # Run Cox sanity check
    coxdir="output/${fam}/all_comorbid/cox/${dim}"
    mkdir -p "$coxdir"
    if [ -f "$coxdir/energy_sanity_check.tsv" ]; then
      echo "  Cox already done, skipping"
    else
      echo "  Running Cox sanity check ..."
      singularity exec docker://koki/desc_onlinexxx_landscaper:20260602 \
        Rscript src/cox_sanity_check_simple.R \
          "$fam" "$dim" \
          "output/${fam}/all_comorbid" \
          "data/tekiyo.csv" \
          "$coxdir" \
        2>&1 | tail -20
      echo "  Cox done"
    fi
  done
done

echo ""
echo "====== Compiling summary ======"
mkdir -p summary
singularity exec docker://koki/desc_onlinexxx_landscaper:20260602 \
  Rscript src/compile_summary.R summary/atlas_summary.tsv 2>&1 | tail -5

echo "====== ALL DONE ======"
