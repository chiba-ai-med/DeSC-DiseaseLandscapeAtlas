# Project Plan: Disease-Family-Specific Landscape Atlas

## Goal

Explore whether disease-family-specific cohort design and feature engineering
can produce energy landscapes where **low-energy basins correspond to severe /
post-complication / high-risk clinical states**.

The key indicator is the Energy Cox sanity check:

- **Pass**: HR per +1 SD Energy < 1 (low energy = high risk)
- **Fail**: HR per +1 SD Energy > 1 (high energy = high risk = frequency landscape)

## Phases

### Phase 0: Skeleton & infrastructure (current)

- Repository structure, config, Snakemake workflows
- 3 pilot disease families defined
- Preprocessing strategies parameterized

### Phase 1: Cohort extraction

For each disease family:

1. Extract patients with at least one claim matching the index ICD-10 codes
2. Identify index date (first diagnosis date)
3. Define baseline window (all claims before index date or within the full
   observation period if temporal splitting is not feasible)
4. Define follow-up window (after index date, for outcome events)
5. Output: `index_cohort.tsv` with patient demographics and observation periods

### Phase 2: Feature matrix construction

For each disease family x preprocessing setting:

1. Collect all disease codes from baseline claims of cohort patients
2. Exclude index disease codes
3. Apply preprocessing filter (all, drop-top-N, complications-only)
4. Build sparse binary COO matrix (patients x disease codes)
5. Output: `feature_matrix.coo.zst`, `row_ids.txt`, `col_ids.txt`

### Phase 3: Dimensionality reduction

1. Run OnlinePCA.jl (exact out-of-core PCA) on each feature matrix
2. Compute for dim = 2 through 7
3. Output: `Scores.mm`, `Eigen_vectors.csv`

### Phase 4: Energy Landscape Analysis

1. Landscaper reads `Scores.mm` (binarized PCA scores, +1/-1)
2. Fits Ising model, identifies basins
3. Output: `Allstates.tsv`, `Basin.tsv`, `SubGraph.tsv`, `E.tsv`

### Phase 5: Energy Cox sanity check

For each setting:

1. Assign each patient a basin and energy value
2. Z-standardize energy within cohort
3. Fit Cox models:
   - `mortality ~ E_z + age + sex`
   - `hospitalization ~ E_z + age + sex`
4. Record HR, 95% CI, p-value
5. Flag `energy_direction_ok = TRUE/FALSE`

### Phase 6: Summary & comparison

1. Compile `atlas_summary.tsv` across all disease x preproc x dim combinations
2. Identify which settings produce the expected energy direction
3. Annotate basin interpretability

### Phase 7: Visualization

1. Landscape plots (Landscaper standard outputs)
2. KM curves by energy tertile
3. Forest plots of energy HR across settings
4. Heatmap: disease x preproc x dim -> energy direction

## Decision Framework

After Phase 5, for each disease family:

- If **all_comorbid** already passes the sanity check -> proceed with full
  analysis using this simple preprocessing
- If only **drop_topN** passes -> determine the optimal N and use that
- If only **complications** passes -> consider whether the code list is too
  restrictive and whether the landscape is clinically meaningful
- If **nothing** passes -> revisit the cohort definition (e.g., narrower index
  disease, different baseline window, add biomarkers)

## Extensibility

- Strategy D (health checkup biomarkers) is deferred to Phase 2+ if code-only
  features do not yield satisfactory landscapes
- Additional disease families can be added by extending
  `config/disease_families.tsv`
