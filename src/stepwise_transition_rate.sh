#!/bin/bash
# Calculate step-by-step monthly basin transition rate for one disease
# Usage: bash stepwise_transition_rate.sh <disease> <output_tsv>
# Input: plot/{disease}/trajectories_dead.tsv
# Output: TSV with from_basin, to_basin, count, rate

DIS=$1
OUT=$2
BASE=$(dirname $(dirname $(dirname "$0")))
TRAJ="${BASE}/plot/${DIS}/trajectories_dead.tsv"
BFILE="${BASE}/output/${DIS}/all_comorbid/landscaper/7/Basin.tsv"

NB=$(wc -w < "$BFILE")

awk -F'\t' -v dis="$DIS" -v nb="$NB" '
NR>1 {
  pid=$1; bid=$5+0; t=$9+0
  if (pid == prev_pid && t == prev_t - 1) {
    key = prev_bid "\t" bid
    trans[key]++
    from_total[prev_bid]++
  }
  prev_pid=pid; prev_bid=bid; prev_t=t
}
END {
  printf "disease\tn_basins\tfrom_basin\tto_basin\tcount\trate\n"
  for (k in trans) {
    split(k, a, "\t")
    printf "%s\t%d\t%d\t%d\t%d\t%.6f\n", dis, nb, a[1], a[2], trans[k], trans[k]/from_total[a[1]]
  }
}' "$TRAJ" > "$OUT"
