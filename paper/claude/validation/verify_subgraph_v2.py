#!/usr/bin/env python3
"""Verify SubGraph state-label mapping — SKIP vs TRUNCATE hypothesis.

Confirmed: SubGraph.tsv rows correspond to non-isolated states in
state-ID order (SKIP), but build_patient_scores.R reads them as
row k = state k (TRUNCATE), causing misassignment after the first gap.

This script quantifies the impact across all 51 diseases.
"""

import csv
import json
from pathlib import Path
from collections import Counter

BASE = Path("/home/koki/dev/DeSC-DiseaseLandscapeAtlas")
OUT_DIR = BASE / "paper/claude/validation"


def read_tsv_ints(path):
    with open(path) as f:
        return [int(line.strip()) for line in f if line.strip()]


def read_tsv_floats(path):
    with open(path) as f:
        return [float(line.strip()) for line in f if line.strip()]


def read_status_network(path):
    matrix = []
    with open(path) as f:
        for line in f:
            matrix.append([int(x) for x in line.strip().split()])
    return matrix


def find_isolated(sn):
    """States (0-indexed) with no 1 in their StatusNetwork row."""
    return [i for i in range(len(sn)) if sum(sn[i]) == 0]


def trace_basin_from_network(sn, energies, basin_minima_1idx):
    """Trace steepest descent to assign basin labels from StatusNetwork."""
    n = len(energies)
    basin_min_set = set(b - 1 for b in basin_minima_1idx)
    basin_label_map = {b - 1: idx + 1 for idx, b in enumerate(basin_minima_1idx)}

    def get_target(i):
        ones = [j for j in range(n) if sn[i][j] == 1]
        return ones[0] if ones else None

    labels = {}
    for start in range(n):
        curr = start
        visited = []
        found = False
        while curr not in set(visited):
            if curr in basin_min_set:
                label = basin_label_map[curr]
                for v in visited:
                    labels[v] = label
                labels[curr] = label
                found = True
                break
            visited.append(curr)
            target = get_target(curr)
            if target is None:
                for v in visited:
                    labels[v] = 0
                found = True
                break
            if energies[target] >= energies[curr]:
                for v in visited:
                    labels[v] = 0
                found = True
                break
            curr = target
        if not found:
            # Cycle without basin minimum — assign by lowest energy in cycle
            # Check if any state in the cycle is adjacent to a basin
            for v in visited:
                labels[v] = 0
    return labels


def analyze_disease(disease):
    ls_dir = BASE / f"output/{disease}/all_comorbid/landscaper/7"
    try:
        subgraph = read_tsv_ints(ls_dir / "SubGraph.tsv")
        energies = read_tsv_floats(ls_dir / "E.tsv")
        freqs = read_tsv_ints(ls_dir / "Freq.tsv")
        basins = read_tsv_ints(ls_dir / "Basin.tsv")
        sn = read_status_network(ls_dir / "StatusNetwork.tsv")
    except Exception as e:
        return {"disease": disease, "error": str(e)}

    n_states = len(energies)
    n_subgraph = len(subgraph)
    n_basins = len(basins)
    isolated = find_isolated(sn)
    n_isolated = len(isolated)
    isolated_set = set(isolated)

    # Non-isolated states in order (SKIP mapping)
    non_isolated = sorted([i for i in range(n_states) if i not in isolated_set])

    # Determine SubGraph format
    if n_subgraph == n_states:
        # Full 128 rows — no skip, SubGraph[k] = state k
        sg_format = "full_128"
        skip_labels = {i: subgraph[i] for i in range(n_states)}
        trunc_labels = {i: subgraph[i] for i in range(n_states)}
    elif len(non_isolated) == n_subgraph:
        sg_format = "skip_isolated"
        # SKIP mapping: SubGraph[k] -> non_isolated[k]
        skip_labels = {}
        for k, state in enumerate(non_isolated):
            skip_labels[state] = subgraph[k]
        for state in isolated:
            skip_labels[state] = 0
        # TRUNCATE mapping (what the code does): SubGraph[k] -> state k
        trunc_labels = {}
        for k in range(n_subgraph):
            trunc_labels[k] = subgraph[k]
        for k in range(n_subgraph, n_states):
            trunc_labels[k] = None  # NA in merge
    else:
        sg_format = f"unknown({len(non_isolated)}ni_vs_{n_subgraph}sg)"
        skip_labels = {i: subgraph[i] if i < n_subgraph else None for i in range(n_states)}
        trunc_labels = skip_labels.copy()

    # Network-derived labels
    network_labels = trace_basin_from_network(sn, energies, basins)

    # Compare SKIP vs network
    skip_vs_network_match = sum(
        1 for i in range(n_states)
        if skip_labels.get(i, -1) == network_labels.get(i, -2)
    )

    # Compare TRUNCATE vs network
    trunc_vs_network = []
    for i in range(n_states):
        t = trunc_labels.get(i)
        n_lab = network_labels.get(i, -1)
        if t is None:
            trunc_vs_network.append("NA_vs_" + str(n_lab))
        elif t == n_lab:
            trunc_vs_network.append("match")
        else:
            trunc_vs_network.append("mismatch")

    # Count misassigned observations
    total_obs = sum(freqs)
    misassigned_states = []
    misassigned_obs = 0
    for i in range(n_states):
        t = trunc_labels.get(i)
        s = skip_labels.get(i, 0)
        if t != s:
            misassigned_states.append(i + 1)  # 1-indexed
            misassigned_obs += freqs[i]

    # Categorize misassignment types
    wrong_basin = 0  # assigned to wrong basin
    false_basin = 0  # isolated assigned to a basin
    lost_basin = 0   # basin state gets NA
    wrong_basin_obs = 0
    false_basin_obs = 0
    lost_basin_obs = 0

    for i in range(n_states):
        correct = skip_labels.get(i, 0)
        assigned = trunc_labels.get(i)
        if correct == assigned:
            continue
        if i in isolated_set and assigned is not None and assigned > 0:
            false_basin += 1
            false_basin_obs += freqs[i]
        elif assigned is None and correct > 0:
            lost_basin += 1
            lost_basin_obs += freqs[i]
        elif assigned is not None and correct > 0 and assigned != correct:
            wrong_basin += 1
            wrong_basin_obs += freqs[i]
        elif assigned is not None and assigned > 0 and correct == 0:
            false_basin += 1
            false_basin_obs += freqs[i]

    return {
        "disease": disease,
        "n_states": n_states,
        "n_subgraph": n_subgraph,
        "n_basins": n_basins,
        "n_isolated": n_isolated,
        "isolated_states": [i + 1 for i in isolated],
        "sg_format": sg_format,
        "skip_network_match": skip_vs_network_match,
        "n_misassigned_states": len(misassigned_states),
        "misassigned_states": misassigned_states,
        "total_obs": total_obs,
        "misassigned_obs": misassigned_obs,
        "misassigned_obs_pct": 100.0 * misassigned_obs / total_obs if total_obs > 0 else 0,
        "wrong_basin_states": wrong_basin,
        "wrong_basin_obs": wrong_basin_obs,
        "false_basin_states": false_basin,
        "false_basin_obs": false_basin_obs,
        "lost_basin_states": lost_basin,
        "lost_basin_obs": lost_basin_obs,
    }


def main():
    diseases_file = BASE / "config/disease_families.tsv"
    diseases = []
    with open(diseases_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            d = row.get('family_id') or row.get('disease_family') or row.get('disease')
            if d and d != 'hf_cancer':
                diseases.append(d)

    results = []
    print(f"{'Disease':<35} {'SubG':>4} {'Isol':>4} {'Bas':>3} {'Format':<16} "
          f"{'MisStates':>9} {'MisObs%':>8} {'WrongBas':>8} {'FalseBas':>8} {'LostBas':>8}")
    print("-" * 135)

    for disease in sorted(diseases):
        r = analyze_disease(disease)
        results.append(r)
        if r.get("error"):
            print(f"{disease:<35} ERROR: {r['error']}")
            continue
        print(f"{disease:<35} {r['n_subgraph']:>4} {r['n_isolated']:>4} {r['n_basins']:>3} {r.get('sg_format','?'):<16} "
              f"{r['n_misassigned_states']:>9} {r['misassigned_obs_pct']:>7.2f}% "
              f"{r['wrong_basin_states']:>8} {r['false_basin_states']:>8} {r['lost_basin_states']:>8}")

    # Save full results
    with open(OUT_DIR / "subgraph_verification_v2.json", 'w') as f:
        json.dump(results, f, indent=2, default=str)

    # Save TSV summary
    with open(OUT_DIR / "subgraph_reconciliation.tsv", 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["disease", "n_subgraph", "n_isolated", "n_basins",
                     "n_misassigned_states", "misassigned_obs_pct",
                     "wrong_basin_states", "wrong_basin_obs",
                     "false_basin_states", "false_basin_obs",
                     "lost_basin_states", "lost_basin_obs",
                     "total_obs", "isolated_states"])
        for r in results:
            if r.get("error"):
                w.writerow([r["disease"], "ERROR"])
                continue
            w.writerow([
                r["disease"], r["n_subgraph"], r["n_isolated"], r["n_basins"],
                r["n_misassigned_states"], f"{r['misassigned_obs_pct']:.2f}",
                r["wrong_basin_states"], r["wrong_basin_obs"],
                r["false_basin_states"], r["false_basin_obs"],
                r["lost_basin_states"], r["lost_basin_obs"],
                r["total_obs"],
                ";".join(str(s) for s in r.get("isolated_states", []))
            ])

    # Summary
    print("\n" + "=" * 80)
    affected = [r for r in results if not r.get("error") and r["n_isolated"] > 0]
    unaffected = [r for r in results if not r.get("error") and r["n_isolated"] == 0]
    print(f"Diseases with isolated states (affected): {len(affected)}")
    print(f"Diseases with 128 SubGraph rows (unaffected): {len(unaffected)}")

    if affected:
        multi_basin_affected = [r for r in affected if r["n_basins"] > 1]
        single_basin_affected = [r for r in affected if r["n_basins"] == 1]
        print(f"  Multi-basin affected: {len(multi_basin_affected)}")
        print(f"  Single-basin affected: {len(single_basin_affected)}")

        total_wrong = sum(r["wrong_basin_states"] for r in affected)
        total_false = sum(r["false_basin_states"] for r in affected)
        total_lost = sum(r["lost_basin_states"] for r in affected)
        print(f"\n  Total wrong-basin assignments: {total_wrong} states")
        print(f"  Total false-basin (isolated→basin): {total_false} states")
        print(f"  Total lost-basin (basin→NA): {total_lost} states")

        # Diseases with wrong_basin > 0 (most serious)
        if multi_basin_affected:
            print(f"\n  Multi-basin diseases with misassignment:")
            for r in sorted(multi_basin_affected, key=lambda x: -x["misassigned_obs_pct"]):
                print(f"    {r['disease']}: {r['n_basins']} basins, "
                      f"{r['n_isolated']} isolated, "
                      f"{r['n_misassigned_states']} mis-states, "
                      f"{r['misassigned_obs_pct']:.2f}% obs, "
                      f"wrong_basin={r['wrong_basin_states']}")


if __name__ == "__main__":
    main()
