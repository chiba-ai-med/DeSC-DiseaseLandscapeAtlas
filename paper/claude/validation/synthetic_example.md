# Synthetic 8-State Example: SubGraph Row-Order vs State-ID Mapping

## Setup

Consider a 3-dimensional Ising model with 2³ = 8 states.
States 3 and 7 are isolated local minima (no lower-energy Hamming-1 neighbor).
State 2 is the basin minimum of Basin A; State 6 is the basin minimum of Basin B.

| State ID | Spin pattern | Energy | Isolated? | True basin |
|----------|-------------|--------|-----------|------------|
| 1        | (−1,−1,−1) | −0.20  | No        | Basin A    |
| 2        | (+1,−1,−1) | −0.80  | No        | Basin A (minimum) |
| 3        | (−1,+1,−1) | −0.15  | **Yes**   | Non-basin  |
| 4        | (+1,+1,−1) | −0.30  | No        | Basin A    |
| 5        | (−1,−1,+1) | −0.25  | No        | Basin B    |
| 6        | (+1,−1,+1) | −0.70  | No        | Basin B (minimum) |
| 7        | (−1,+1,+1) | −0.10  | **Yes**   | Non-basin  |
| 8        | (+1,+1,+1) | −0.35  | No        | Basin B    |

## Landscaper Output: SubGraph.tsv (6 rows, states 3 and 7 excluded)

```
1     ← state 1 → Basin A
1     ← state 2 → Basin A
1     ← state 4 → Basin A
2     ← state 5 → Basin B
2     ← state 6 → Basin B
2     ← state 8 → Basin B
```

## What the code does (TRUNCATE: row k = state k)

```r
subgraph[, pattern_id := .I]  # pattern_id = 1,2,3,4,5,6
```

| Patient state | pattern_id | SubGraph row matched | Assigned basin_id | Correct basin | Match? |
|:------------:|:----------:|:-------------------:|:-----------------:|:------------:|:------:|
| 1            | 1          | Row 1 (=state 1)    | 1 (Basin A)       | Basin A      | ✓      |
| 2            | 2          | Row 2 (=state 2)    | 1 (Basin A)       | Basin A      | ✓      |
| 3            | 3          | Row 3 (=state **4**)| 1 (Basin A)       | **Non-basin**| **✗**  |
| 4            | 4          | Row 4 (=state **5**)| 2 (Basin B)       | **Basin A**  | **✗**  |
| 5            | 5          | Row 5 (=state **6**)| 2 (Basin B)       | Basin B      | ✓ (luck) |
| 6            | 6          | Row 6 (=state **8**)| 2 (Basin B)       | Basin B      | ✓ (luck) |
| 7            | 7          | — (no match)        | NA                | **Non-basin**| ✓ (luck) |
| 8            | 8          | — (no match)        | NA                | **Basin B**  | **✗**  |

**Result**: 4/8 states misassigned (50%). State 3 (isolated) gets basin A; state 4 gets wrong basin; state 8 loses its basin.

## What the correct mapping would be (SKIP: row k = non-isolated state k)

SubGraph rows map to non-isolated states in state-ID order:
Row 1→State 1, Row 2→State 2, Row 3→State 4, Row 4→State 5, Row 5→State 6, Row 6→State 8

| Patient state | pattern_id | Correct SubGraph row | Correct basin_id | Match? |
|:------------:|:----------:|:-------------------:|:----------------:|:------:|
| 1            | 1          | Row 1               | 1 (Basin A)      | ✓      |
| 2            | 2          | Row 2               | 1 (Basin A)      | ✓      |
| 3            | 3          | — (isolated)        | 0 (Non-basin)    | ✓      |
| 4            | 4          | Row 3               | 1 (Basin A)      | ✓      |
| 5            | 5          | Row 4               | 2 (Basin B)      | ✓      |
| 6            | 6          | Row 5               | 2 (Basin B)      | ✓      |
| 7            | 7          | — (isolated)        | 0 (Non-basin)    | ✓      |
| 8            | 8          | Row 6               | 2 (Basin B)      | ✓      |

**Result**: 8/8 states correctly assigned.

## Implication for build_patient_scores.R

Line 36: `subgraph[, pattern_id := .I]` assigns sequential row numbers as pattern_id.
Line 69: `merge(scores, subgraph, by = "pattern_id", all.x = TRUE)` joins by pattern_id.

Since patient pattern_ids are 1–128 (from Allstates.tsv) but SubGraph pattern_ids are 1–N
(N < 128), the join is correct only for states before the first isolated-state gap.
After the gap, every state is assigned the wrong SubGraph row's label.
