# Paper Figure Candidates

## Figure Structure

- **Fig 1**: Schematic diagram of ELA pipeline (to be drawn manually)
- **Fig 2**: Atlas heatmap — PASS/FAIL across 51 disease families
- **Fig 3**: Applicability scatter — death rate vs observation density (PASS/FAIL colored)
- **Fig 4**: Representative disease comparison (diabetes, dementia vs heart_failure, sepsis)
  - 4a: Energy trajectory before death
  - 4b: KM survival by energy tertile (most informative PASS/FAIL contrast)
  - 4c: State occupancy heatmap (top 10 states)
- **Fig 5**: Trajectory-to-death summary (all 51 diseases)
  - 5a: Non-basin fraction at death (bar chart) — median 4.8%, max 13.2%
  - 5b: Mean entropy trajectory PASS vs FAIL
  - 5c: Entropy drop vs death rate scatter
- **Fig 6**: Per-basin trajectory (representative diseases)
  - 6a: Dementia energy trajectory per basin (2 basins, clear separation)
  - 6b: Diabetes energy trajectory per basin (4 basins, convergence at death)
  - 6c: Dementia Basin + Non-basin Sankey (Non-basin ~0%)
  - 6d: Diabetes Basin + Non-basin Sankey (Non-basin ~0%)

## Supplementary Figures

- `FigS_dementia_basin_sankey.png` — Basin-to-basin transitions only (dementia)
- `FigS_diabetes_basin_sankey.png` — Basin-to-basin transitions only (diabetes)
- `FigS_heart_failure_sankey_nonbasin.png` — FAIL example: heart_failure
- `FigS_sepsis_sankey_nonbasin.png` — FAIL example: sepsis
- `FigS_entropy_all_51.png` — All 51 diseases entropy trajectories overlaid
- `FigS_n_states_at_death.png` — Number of distinct states at death (bar chart)
- Per-disease plots available in `plot/{disease}/` for all 51 diseases
