#!/bin/bash
## Build window lookup for all diseases at once.
## 1. Collect all old_row_ids from dead patients across all diseases
## 2. Single pass through row_id_number_small.txt to extract matches
## 3. Output: data/rowid_window_lookup.tsv (row_id \t window)
##
## Usage: bash src/build_window_lookup.sh

set -e

echo "=== Step 1: Collect dead patient old_row_ids ==="
> /tmp/dead_old_row_ids.txt
for d in $(awk -F'\t' 'NR>1{print $1}' config/disease_families.tsv); do
  pt="output/$d/all_comorbid/patient/patient_scores.tsv"
  pmap="output/$d/all_comorbid/cohort/person_map.tsv"
  [ -f "$pt" ] || continue
  [ -f "$pmap" ] || continue

  # Get kojin_ids of dead patients
  awk -F'\t' 'NR>1 && $15=="dead"{print $1}' "$pt" > /tmp/dead_kojin_$d.txt
  n_dead=$(wc -l < /tmp/dead_kojin_$d.txt)

  # Get their old_row_ids from person_map
  awk -F'\t' 'NR==FNR{dead[$1]; next} $3 in dead{print $2}' \
    /tmp/dead_kojin_$d.txt "$pmap" >> /tmp/dead_old_row_ids.txt

  echo "  $d: $n_dead dead patients"
  rm /tmp/dead_kojin_$d.txt
done

sort -un /tmp/dead_old_row_ids.txt > /tmp/dead_old_row_ids_uniq.txt
n_total=$(wc -l < /tmp/dead_old_row_ids_uniq.txt)
echo "  Total unique old_row_ids: $n_total"

echo "=== Step 2: Single pass through row_id_number ==="
awk 'NR==FNR{ids[$1]; next} $2 in ids{split($1,a,"_"); print $2 "\t" a[1]}' \
  /tmp/dead_old_row_ids_uniq.txt data/row_id_number_small.txt \
  > data/rowid_window_lookup.tsv

n_matched=$(wc -l < data/rowid_window_lookup.tsv)
echo "  Matched: $n_matched rows"
echo "  Written: data/rowid_window_lookup.tsv"

rm /tmp/dead_old_row_ids.txt /tmp/dead_old_row_ids_uniq.txt
echo "Done."
