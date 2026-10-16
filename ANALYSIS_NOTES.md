# Analysis Notes

## What was done

### 1. Trajectory-to-death analysis (all 51 diseases)
**Script**: `src/trajectory_to_death.R`
**Output**: `plot/{disease}/trajectories_dead.tsv`, `trajectory_heatmap.png`, etc.

For each deceased patient, track their state (Allstates pattern), energy, and basin
membership at every observation window, indexed by periods-before-death (up to 20
periods). Produces per-disease heatmaps of state occupancy over time before death.

### 2. Applicability analysis: why PASS vs FAIL?
**Script**: `src/applicability_analysis.R`
**Output**: `summary/applicability/disease_metrics.tsv`, scatter plots

Quantified what distinguishes 43 PASS diseases from 8 FAIL diseases.
- Logistic regression: `PASS ~ death_rate + log(obs_per_patient)`
  - death_rate: p=0.028
  - log(obs_per_patient): p=0.007
  - AIC = 34.3
- FAIL diseases tend to have lower observation density and/or acute presentation
  (sepsis, delirium, acute_renal_failure, tuberculosis)
- Interpretation: PASS = chronic progressive diseases where basin = severity stage.
  FAIL = acute/episodic diseases where instability (high energy) = risk.

### 3. Representative disease comparison
**Scripts**: `src/figure_representative_comparison.R`, `src/figure_km_comparison.R`
**Output**: `summary/applicability/panel_*.png`

Compared 2 PASS (diabetes HR=0.97, dementia HR=0.99) vs 2 FAIL (heart_failure
HR=1.04, sepsis HR=1.06):
- Energy trajectory: PASS shows declining energy before death (patients sink into
  low-energy basins). FAIL shows flat or increasing energy.
- KM by energy tertile: PASS shows blue (low energy) below red (high energy) = low
  energy worse survival. FAIL shows reversed or no separation.
- State heatmap: PASS shows concentration into fewer states. FAIL shows dispersed.

### 4. Per-basin trajectory analysis
**Script**: `src/trajectory_per_basin.R`
**Output**: `plot/{disease}/trajectory_energy_per_basin.png`,
`trajectory_basin_composition.png`, `trajectory_basin_sankey.png`

For each disease with 2+ basins:
- Energy trajectory faceted by basin: shows that different basins have distinct
  energy levels and trajectories.
- Basin composition: stacked area chart showing fraction of patients in each basin
  over time before death. Very stable — patients stay in their basin.
- Basin transition Sankey: ~94% stay in same basin per period, ~5% transition.

### 5. State entropy + Non-basin Sankey
**Script**: `src/trajectory_entropy_sankey.R`
**Output**: `plot/{disease}/trajectory_entropy.png`, `trajectory_n_states.png`,
`trajectory_sankey_with_nonbasin.png`

- Shannon entropy (normalized 0-1) across 128 states at each time point before
  death. Measures how concentrated patients are across states.
- Basin + Non-basin Sankey: includes Non-basin (basin_id=0) as a node category.

### 6. Entropy summary across all 51 diseases
**Script**: `src/plot_entropy_summary.R`
**Output**: `summary/applicability/entropy_*.png`, `barplot_*.png`,
`scatter_entropy_*.png`, `entropy_summary.tsv`

Key findings:
- **Non-basin fraction at death**: median 4.8%, range 0-13.2%. Patients are
  overwhelmingly in basins at time of death. Lowest: sepsis/diabetes/dementia
  (~0%). Highest: malnutrition (13.2%), tuberculosis (9.1%).
- **Entropy**: all diseases 0.93-0.98 (high). Slight decrease before death.
  No major PASS/FAIL difference. Cohort-level entropy stays high because 128
  states are well-populated in large cohorts.
- **Entropy drop**: PASS diseases show small positive drop (patients concentrate
  slightly before death). FAIL diseases sepsis/tuberculosis show large negative
  drop (patients disperse before death).

## Key results summary

1. **ELA works for 43/51 diseases** (84%): low energy = high mortality risk (HR<1).
2. **PASS/FAIL predictable** from death_rate + observation density (both significant).
3. **Patients stay in basins**: 87-100% at death (median 95.2%). Basin transitions
   are rare (~5% per period).
4. **Basin = stable comorbidity pattern**: patients don't leave basins before death.
   This validates basins as disease severity stages.
5. **FAIL diseases** (sepsis, delirium, etc.) are acute/episodic — instability
   rather than low energy predicts risk.

## Snakemake

All new scripts are registered in `workflow/figures.smk`:
- `trajectory_per_basin`
- `trajectory_entropy_sankey`
- `applicability_analysis`
- `figure_representative_comparison`
- `figure_km_comparison`
- `plot_entropy_summary`

## Source files added/modified

- `src/trajectory_per_basin.R` — new
- `src/trajectory_entropy_sankey.R` — new
- `src/plot_entropy_summary.R` — new
- `src/applicability_analysis.R` — new
- `src/figure_representative_comparison.R` — new
- `src/figure_km_comparison.R` — new
- `workflow/figures.smk` — updated (6 rules added)
- `plot/Figures/` — paper figure candidates
