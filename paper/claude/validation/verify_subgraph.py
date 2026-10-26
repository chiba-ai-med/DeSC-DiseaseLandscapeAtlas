#!/usr/bin/env python3
"""Verify SubGraph state-label mapping for all 51 diseases.

For each disease:
1. Read Allstates.tsv (128 states), E.tsv (128 energies), Basin.tsv (local minima)
2. Reconstruct steepest-descent basin assignment from StatusNetwork.tsv
3. Compare with SubGraph.tsv (which may have <128 rows)
4. Identify the mapping between SubGraph rows and state IDs
"""

import os
import sys
import csv
from pathlib import Path
import json

BASE = Path("/home/koki/dev/DeSC-DiseaseLandscapeAtlas")
OUT_DIR = BASE / "paper/claude/validation"


def read_allstates(disease):
    """Read 128 spin patterns."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/Allstates.tsv"
    states = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            states.append(tuple(int(x) for x in parts))
    return states


def read_energies(disease):
    """Read 128 energies."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/E.tsv"
    with open(path) as f:
        return [float(line.strip()) for line in f]


def read_frequencies(disease):
    """Read 128 frequencies."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/Freq.tsv"
    with open(path) as f:
        return [int(line.strip()) for line in f]


def read_basins(disease):
    """Read basin local minima state IDs (1-indexed)."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/Basin.tsv"
    with open(path) as f:
        return [int(line.strip()) for line in f]


def read_subgraph(disease):
    """Read SubGraph labels."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/SubGraph.tsv"
    with open(path) as f:
        return [int(line.strip()) for line in f]


def read_status_network(disease):
    """Read 128x128 adjacency matrix from StatusNetwork.tsv."""
    path = BASE / f"output/{disease}/all_comorbid/landscaper/7/StatusNetwork.tsv"
    matrix = []
    with open(path) as f:
        for line in f:
            row = [int(x) for x in line.strip().split()]
            matrix.append(row)
    return matrix


def compute_steepest_descent_basins(energies, status_net):
    """
    From StatusNetwork adjacency, compute basin assignment for each state.
    StatusNetwork[i] has a 1 at position j if state i connects to state j
    in steepest descent. Follow the chain until reaching a local minimum
    (a state that points to itself or has no outgoing edge).
    """
    n = len(energies)
    # For each state, find the neighbor with lower energy
    # StatusNetwork rows encode the graph structure
    # Each row has 1s at positions of connected states

    # Find steepest-descent neighbor for each state
    next_state = list(range(n))  # default: self
    for i in range(n):
        neighbors = [j for j in range(n) if status_net[i][j] == 1]
        if not neighbors:
            next_state[i] = i  # isolated or local minimum
            continue
        # Include self in comparison
        candidates = neighbors + [i]
        best = min(candidates, key=lambda x: energies[x])
        next_state[i] = best

    # Follow chains to local minima
    basin_of = [None] * n
    for i in range(n):
        visited = []
        curr = i
        while basin_of[curr] is None:
            if curr in visited:
                # cycle - assign to lowest energy in cycle
                cycle_min = min(visited, key=lambda x: energies[x])
                for v in visited:
                    basin_of[v] = cycle_min
                break
            visited.append(curr)
            nxt = next_state[curr]
            if nxt == curr:
                # local minimum
                for v in visited:
                    basin_of[v] = curr
                break
            curr = nxt
        else:
            # curr already assigned
            for v in visited:
                basin_of[v] = basin_of[curr]

    return basin_of


def find_local_minima(energies, status_net):
    """Find states that are local minima (no lower-energy neighbor)."""
    n = len(energies)
    minima = []
    for i in range(n):
        neighbors = [j for j in range(n) if status_net[i][j] == 1]
        is_minimum = all(energies[j] >= energies[i] for j in neighbors)
        if is_minimum:
            minima.append(i)
    return minima


def analyze_disease(disease):
    """Full analysis for one disease."""
    try:
        allstates = read_allstates(disease)
        energies = read_energies(disease)
        freqs = read_frequencies(disease)
        basins = read_basins(disease)  # 1-indexed state IDs
        subgraph = read_subgraph(disease)
        status_net = read_status_network(disease)
    except Exception as e:
        return {"disease": disease, "error": str(e)}

    n_states = len(allstates)
    n_subgraph = len(subgraph)
    n_basins = len(basins)

    # Basin state IDs are 1-indexed in Basin.tsv
    basin_states_0idx = [b - 1 for b in basins]

    # Find local minima from energy landscape
    local_minima = find_local_minima(energies, status_net)

    # Check if Basin.tsv minima match actual local minima
    basin_match_minima = set(basin_states_0idx) == set(local_minima)

    # Compute expected basin assignment from steepest descent
    sd_basin = compute_steepest_descent_basins(energies, status_net)

    # Map basin assignments: sd_basin gives the state_id of the minimum
    # Convert to basin labels (1-indexed matching Basin.tsv order)
    basin_label_map = {basin_states_0idx[i]: i + 1 for i in range(n_basins)}

    expected_labels = []
    for i in range(n_states):
        if sd_basin[i] in basin_label_map:
            expected_labels.append(basin_label_map[sd_basin[i]])
        else:
            expected_labels.append(0)  # non-basin

    # Determine what SubGraph actually represents
    # Hypothesis 1: SubGraph has one row per state that is NOT a singleton
    # (singleton = isolated component in the status network = freq too low)
    # Hypothesis 2: SubGraph is indexed by state_id for non-removed states

    # Find which states might be singletons (no connections)
    connected_states = []
    for i in range(n_states):
        neighbors = [j for j in range(n_states) if status_net[i][j] == 1]
        if neighbors:
            connected_states.append(i)

    # Check if n_subgraph == n_connected
    n_connected = len(connected_states)

    # Also check by frequency - singletons might be zero-frequency states
    zero_freq = [i for i in range(n_states) if freqs[i] == 0]
    nonzero_freq = [i for i in range(n_states) if freqs[i] > 0]

    # Test hypotheses
    # H1: SubGraph rows correspond to connected states (in order)
    h1_match = n_subgraph == n_connected
    # H2: SubGraph rows correspond to non-zero-frequency states
    h2_match = n_subgraph == len(nonzero_freq)

    # For TB: 128 states, 118 SubGraph rows, so 10 states removed
    removed_states = n_states - n_subgraph
    removed_by_connection = n_states - n_connected
    removed_by_freq = len(zero_freq)

    # Build the actual mapping if H1 or H2 holds
    if h1_match:
        # SubGraph[k] corresponds to connected_states[k]
        actual_labels = [0] * n_states
        for k, state_id in enumerate(connected_states):
            actual_labels[state_id] = subgraph[k]
        mapping_type = "connected_states"
        mapped_states = connected_states
    elif h2_match:
        actual_labels = [0] * n_states
        for k, state_id in enumerate(nonzero_freq):
            actual_labels[state_id] = subgraph[k]
        mapping_type = "nonzero_freq"
        mapped_states = nonzero_freq
    else:
        # Try: first n_subgraph states
        actual_labels = [0] * n_states
        for k in range(n_subgraph):
            actual_labels[k] = subgraph[k]
        mapping_type = "first_n (fallback)"
        mapped_states = list(range(n_subgraph))

    # Compare expected vs actual
    mismatches = []
    for i in range(n_states):
        if expected_labels[i] != actual_labels[i]:
            mismatches.append({
                "state_id": i + 1,  # 1-indexed
                "energy": energies[i],
                "freq": freqs[i],
                "expected_label": expected_labels[i],
                "actual_label": actual_labels[i],
                "is_basin_minimum": i in basin_states_0idx,
            })

    # Compute impact: what fraction of observations are in mismatched states
    total_obs = sum(freqs)
    mismatch_obs = sum(freqs[m["state_id"] - 1] for m in mismatches)

    return {
        "disease": disease,
        "n_states": n_states,
        "n_subgraph": n_subgraph,
        "n_basins": n_basins,
        "basin_states_1idx": basins,
        "n_local_minima": len(local_minima),
        "local_minima_1idx": [m + 1 for m in local_minima],
        "basin_match_minima": basin_match_minima,
        "n_connected": n_connected,
        "n_zero_freq": len(zero_freq),
        "zero_freq_states_1idx": [z + 1 for z in zero_freq],
        "removed_states": removed_states,
        "mapping_type": mapping_type,
        "n_mismatches": len(mismatches),
        "mismatch_obs_frac": mismatch_obs / total_obs if total_obs > 0 else 0,
        "total_obs": total_obs,
        "mismatch_obs": mismatch_obs,
        "mismatches": mismatches[:20],  # first 20 for brevity
    }


def main():
    # Get disease list
    diseases_file = BASE / "config/disease_families.tsv"
    diseases = []
    with open(diseases_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            d = row.get('disease_family') or row.get('disease')
            if d and d != 'hf_cancer':
                diseases.append(d)

    # First, detailed analysis of TB
    print("=" * 60)
    print("DETAILED ANALYSIS: tuberculosis")
    print("=" * 60)
    tb = analyze_disease("tuberculosis")
    for k, v in tb.items():
        if k != "mismatches":
            print(f"  {k}: {v}")
    if tb.get("mismatches"):
        print(f"\n  First mismatches:")
        for m in tb["mismatches"][:10]:
            print(f"    state {m['state_id']}: expected={m['expected_label']}, "
                  f"actual={m['actual_label']}, freq={m['freq']}, E={m['energy']:.4f}")

    # Then all diseases
    print("\n" + "=" * 60)
    print("ALL 51 DISEASES SUMMARY")
    print("=" * 60)

    all_results = []
    for disease in sorted(diseases):
        result = analyze_disease(disease)
        all_results.append(result)
        status = "OK" if result.get("n_mismatches", -1) == 0 else f"MISMATCH({result.get('n_mismatches', '?')})"
        if result.get("error"):
            status = f"ERROR: {result['error']}"
        print(f"  {disease:40s} SubGraph={result.get('n_subgraph','?'):>3} "
              f"basins={result.get('n_basins','?')} "
              f"mapping={result.get('mapping_type','?'):20s} "
              f"{status}")

    # Save full results
    out_path = OUT_DIR / "subgraph_verification.json"
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nFull results saved to {out_path}")

    # Summary stats
    n_ok = sum(1 for r in all_results if r.get("n_mismatches", -1) == 0)
    n_mismatch = sum(1 for r in all_results if r.get("n_mismatches", 0) > 0)
    n_error = sum(1 for r in all_results if r.get("error"))
    print(f"\nSummary: {n_ok} OK, {n_mismatch} with mismatches, {n_error} errors")

    # Diseases with mismatches
    if n_mismatch > 0:
        print("\nDiseases with mismatches:")
        for r in all_results:
            if r.get("n_mismatches", 0) > 0:
                print(f"  {r['disease']}: {r['n_mismatches']} states, "
                      f"mapping={r['mapping_type']}, "
                      f"impact={r['mismatch_obs_frac']:.4%} of observations")


if __name__ == "__main__":
    main()
