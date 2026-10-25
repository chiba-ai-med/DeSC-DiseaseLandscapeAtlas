# Paper Figures — Disease Landscape Atlas (Patterns)

Updated: 2026-10-23

## Figure Structure (9 figures + Table 1)

| Fig/Table | Disease/Topic | Pattern | Key message |
|-----------|--------------|---------|-------------|
| 1 | Workflow | — | ICD-10 → PCA → binarize → 128 states → basins → Cox |
| 2 | PCA (all 51) | — | PC1=age, PC7=sex; 7 PCs capture median 12.5% of total variance |
| 3 | Atlas overview | — | HR < 1 in 43/51 (mortality), 40/51 (hospitalization), per-basin majority |
| 4 | Jaccard clustering | — | Basin code similarity clusters diseases into clinically valid groups |
| 5 | Trajectory | — | Basins trap patients (96% stay); 2+ basin show 50% state transitions |
| Table 1 | 6-pattern classification | All | 51 diseases classified by 4 axes into 6 behavioral patterns |
| 6 | Thyroid | Convergent | 4 basins; CKD complications + textbook hypothyroid symptoms; moderate asymmetry |
| 7 | Dementia | Symmetric | 2 basins; autonomic vs vascular; mortality-hospitalization divergence |
| 8 | Anemia | Progressive | Transition trio comparison; extreme asymmetry (45:1 one-way CKD progression) |
| 9 | TB | Non-predictive | Sankey + PC sign pattern; ELA limitations for low-density acute diseases |

Homogeneous (13 diseases) and Structured non-predictive (2: HF, Sepsis) are text-only.
Diabetes → supplementary (overlaps Thyroid structure). Schizophrenia dropped. Sepsis → supplementary.

## Narrative flow

```
Fig 1-2  Method (workflow + PCA)
Fig 3    Atlas overview (全体俯瞰: energy-mortality scatter)
Fig 4    Jaccard clustering (外的妥当性: basin code類似度で疾患クラスタ再現)
Fig 5    Trajectory (Basin動態)
Table 1  6-pattern classification (4軸 × 51疾患)
  → Homogeneous (13): text only, 1 basin, Cox predictive
  → Convergent:
Fig 6      Thyroid (代表例: bar chart + KM + discon + sankey)
  → Symmetric:
Fig 7      Dementia (代表例: paired KM mortality vs hospitalization)
  → Progressive:
Fig 8      Anemia (代表例: transition trio comparison)
  → Non-predictive:
Fig 9      TB (代表例: Sankey + PC sign pattern)
  → Structured non-predictive (HF, Sepsis): text only, supplementary
```

---

## Figure Titles and Captions

### Fig 1: Overview of the DeSC Disease Landscape Atlas analysis workflow

**Caption:** Overview of the analysis workflow. (a) From the DeSC claims
database, 51 disease-family cohorts are defined by ICD-10 index codes. For each
cohort, a binary comorbidity matrix is constructed (rows: patient-observation
windows, columns: ICD-10 codes excluding the index disease). PCA reduces the
matrix to 7 dimensions, scores are binarized to +/-1, yielding 2^7 = 128 discrete
states. Energy landscape analysis assigns each state an energy value and groups
states into basins (2-5 per disease). (b) Within each basin, Cox proportional
hazards models test whether lower energy is associated with higher mortality or
hospitalization risk (HR < 1). The per-basin majority criterion is satisfied in
43/51 families for mortality and 40/51 for hospitalization.

- `Fig01_workflow.png` — ChatGPT draft (to be traced in Illustrator)

### Fig 2: PCA captures age and sex axes across disease families

**Caption:** (a) Diabetes cohort (n = 4,514,492) projected onto PC1 vs PC2,
colored by age. PC1 correlates with age across most diseases (median r = 0.40).
(b) Urinary tract infection cohort projected onto PC1 vs PC7, colored by sex.
PC7 separates female and male patients (r = -0.45). (c) Cumulative variance
explained by the first 7 PCs for all 51 disease families (grey) and their median
(red). Seven PCs capture a median of 12.5% of total variance; approximately 300
dimensions would be needed for 80%, but ELA requires 2^dim states, making 7 the
practical maximum. (d) Pearson correlation of each PC with age (left) and sex
(right) across 51 disease families, sorted by PC1-age correlation. PC1
universally captures an age gradient; sex associations are distributed across
different PCs per disease.

- `Fig02_pca.png` — Combined figure (assembled externally)
- `Fig02a_diabetes_PC1_PC2_age.png` — (a) Diabetes PC1 vs PC2, age-colored
- `Fig02b_uti_PC1_PC7_sex.png` — (b) UTI PC1 vs PC7, sex-colored
- `Fig02c_cumulative_variance.png` — (c) True cumulative variance explained
- `Fig02d_pc_correlation_vertical.png` — (d) PC-Age/Sex correlation heatmap

### Fig 3: Atlas of 51 disease families by observation density, death rate, and energy-mortality association

**Caption:** Each point represents one disease family. The x-axis shows the
mean number of claims records per patient (observation density), and the y-axis
shows the cohort death rate during follow-up. Point color indicates whether
energy is associated with mortality in the expected direction (HR < 1) in a
majority of basins (blue) or not (red). Point size reflects the magnitude of
the association (|log(median HR)| across basins). The per-basin majority
criterion is satisfied in 43 of 51 disease families. Diseases with low
observation density and high terminality (upper left) are more likely to show
HR >= 1, consistent with insufficient comorbidity information for landscape
construction.

- `Fig03_atlas_overview.png` — Mortality (HR < 1 in 43/51, per-basin majority)
- `FigS_atlas_hospitalization.png` — Supplementary: Hospitalization (HR < 1 in 40/51)

### Fig 4: Disease clustering by basin code similarity recovers clinically valid disease groups

**Caption:** (a) Hierarchical clustering of 50 disease families by basin code
similarity. For each pair of disease families, the maximum Jaccard similarity
between any two basins (using each basin's top-30 enriched ICD-10 codes) is
computed. Average-linkage clustering with optimal leaf ordering reveals two
major groups: an Acute/Frailty cluster (7 diseases: acute renal failure,
sepsis, delirium, aspiration pneumonia, pneumonia, malnutrition, hip fracture)
and a broad Chronic/Lifestyle cluster (43 diseases). Within the Chronic/
Lifestyle cluster, clinically coherent sub-clusters emerge without supervision:
neurodegenerative diseases (Dementia, Epilepsy, Parkinson), prostate diseases
(BPH, Prostate Cancer), pulmonary fibrosis diseases (ILD, Lung Cancer),
musculoskeletal diseases (Osteoarthritis, Osteoporosis), and others. These
groupings, derived solely from ELA basin codes, recapitulate known disease
relationships, providing external validation that ELA basins capture
clinically meaningful comorbidity structure.

- `Fig04_jaccard_clustering.png` — Dendrogram + heatmap (2 large blocks annotated)
- `FigS_jaccard_allblocks.png` — Supplementary: same with 12 sub-clusters annotated

### Fig 5: Patient trajectories reveal basins as traps with basin-count-dependent transition rates

**Caption:** (a) Distribution of dead patients' trajectory patterns by number
of basins. In 1-basin diseases (17 families), 86% of patients remain in the
same basin throughout follow-up. In 2+ basin diseases (34 families), roughly
half undergo state transitions before death. Patients staying entirely outside
basins (Stable Non-Basin) account for less than 2% across all groups.
(b) Among patients starting in a basin, 96% remain in the same basin in
1-basin diseases; this drops to ~68% in 3-4 basin diseases, where 31% move
to a different basin. Exit from basins to non-basin states is rare (1-4%),
consistent with basins acting as energy minima that trap patient trajectories.
(c) Among patients starting outside a basin, 60-70% are absorbed into a basin
before death, while 30-40% remain outside. The absorption rate is slightly
higher in diseases with more basins. In some diseases (e.g. dementia, diabetes,
sepsis), all dead patients start in a basin, so these diseases are absent from
panel (c) (n = 24 instead of 26 for the 2-basin group, etc.). These patterns
are consistent with the energy landscape interpretation: basins are stable
states that attract and retain patients, with inter-basin transitions governed
by barrier heights.

- `Fig05_trajectory.png` — Combined figure (assembled externally)
- `Fig05a_trajectory_overview.png` — (a) Stable Basin / Stable Non-Basin / Moved
- `Fig05a_trajectory_overview.tsv` — Summary data
- `Fig05b_trajectory_basin_start.png` — (b) Basin start → same / different / non-basin
- `Fig05b_trajectory_basin_start.tsv` — Summary data
- `Fig05b_trajectory_basin_start_detail.tsv` — Per-disease detail
- `Fig05c_trajectory_nonbasin_start.png` — (c) Non-Basin start → basin / stay
- `Fig05c_trajectory_nonbasin_start.tsv` — Summary data
- `Fig05c_trajectory_nonbasin_start_detail.tsv` — Per-disease detail

### Fig 6: Thyroid disorder energy landscape rediscovers textbook hypothyroid complications and reveals basin dynamics

**Caption:** (a) Top 10 comorbidity codes enriched in Basin 3 (n = 984,461)
and Basin 4 (n = 139,993) of the thyroid disorder cohort, ranked by log-odds
ratio. Bar color indicates the frequency of each code within the basin. Basin 3
is characterized by chronic kidney disease Stage 5 (N185, 1.8%), hyperkalemia
(E875, 4.0%), renal failure (N19, 7.7%), secondary hyperparathyroidism (E211,
2.2%), and chronic disease anemia (D638, 6.4%), consistent with the known
pathway from hypothyroidism to CKD progression. Constipation (49.2%) and
xerosis cutis (22.7%), classic hypothyroid symptoms, are also highest in
Basin 3 among all four basins. Basin 4 is characterized almost entirely by
dental codes (chronic periodontitis 80.0%, dental caries 34.9%), representing
young patients with active dental care. (b) Kaplan-Meier survival curves for
all four basins (log-rank p < 0.0001). The mortality gap between Basin 2
(n = 120,160; mean age 74.7, death rate 50.4%) and Basin 4 (mean age 61.8,
death rate 30.0%) is 20.4 percentage points, the largest across all 51
disease families. (c) Disconnectivity graph showing the hierarchical basin
structure. Basin 1 and Basin 3 merge at the lowest barrier (E = -0.68),
indicating that these CKD-enriched basins are energetically adjacent.
Basin 2 and Basin 4 join at higher barriers, reflecting greater separation.
(d) Temporal flow of basin composition among deceased patients, shown as
months before death. Basin proportions remain largely stable over the 20-month
window, consistent with basins acting as traps that retain patients.

- `Fig06_thyroid.png` — Combined figure (assembled externally)
- `Fig06a_thyroid_basin3_vs_4.png` — (a) Basin 3 vs 4 top-10 log-odds ratio bar chart
- `Fig06b_thyroid_km_inter_basin.png` — (b) Inter-basin KM (survminer, p < 0.0001)
- `Fig06c_thyroid_disorder_discon.png` — (c) Disconnectivity graph (4 basins)
- `Fig06d_thyroid_disorder_sankey.png` — (d) Temporal basin flow (months before death)

### Fig 7: Dementia — same mortality but divergent hospitalization between autonomic dysfunction and vascular risk basins (Symmetric)

**Caption:** (a) Kaplan-Meier curves for mortality (left) and hospitalization
(right) in the dementia cohort (n = 929,951), stratified by basin. The two
basins show nearly identical mortality (Basin 1: 49.2%, Basin 2: 53.8%;
curves overlap visually despite log-rank p < 0.0001 driven by large sample
size), but strikingly different hospitalization patterns (Basin 1: 77.6%,
Basin 2: 57.2%). Basin 1 (mean age 83.0) is enriched for constipation (79%),
xerosis cutis (39%), insomnia (47%), and heart failure (39%), consistent with
autonomic dysfunction and advanced frailty. Basin 2 (mean age 80.8) is
enriched for hypertension (78%), diabetes (48%), and hyperlipidemia (48%),
consistent with a vascular risk factor burden. The mortality-hospitalization
divergence — where basins with similar death rates differ in healthcare
utilization — suggests that ELA captures care-intensity dimensions orthogonal
to terminal prognosis. The Symmetric transition pattern (Fig 7a) confirms
that these two basins exchange patients at equal rates (5.7% vs 5.9%/month),
consistent with subtype coexistence rather than disease progression.

- `Fig07_dementia.png` — Combined figure (assembled externally)
- `Fig07a_dementia_paired_km.png` — (a) Paired mortality vs hospitalization KM

### Fig 8: Transition dynamics across three ELA behavioral patterns (Anemia as Progressive representative)

**Caption:** Step-by-step monthly basin transition matrices for three
representative diseases, ordered by increasing transition asymmetry.
(a) Dementia (Symmetric pattern): Basin 1→2 = 5.7%, Basin 2→1 = 5.9%,
near-perfect symmetry indicating subtype-like oscillation between basins.
(b) Thyroid disorder (Convergent pattern): 4×4 matrix with moderate
asymmetry; Basin 2→1 (10.2%) exceeds Basin 1→2 (2.3%), indicating net
flow toward Basin 1. Basin 3 is the most stable (91.3% retention).
(c) Anemia (Progressive pattern): Basin 3→1 = 22.6% vs Basin 1→3 = 0.5%,
a 45:1 asymmetry ratio. Basin 1 (CKD/dialysis, E = -0.990) absorbs
patients from Basin 3 (dental/mild, E = +0.146) in a one-way progression
consistent with CKD-anemia disease course.

- `Fig08_anemia.png` — Combined figure (assembled externally)
- `Fig08a_transition_dementia.png` — (a) Dementia 2×2 transition matrix
- `Fig08b_transition_thyroid.png` — (b) Thyroid 4×4 transition matrix
- `Fig08c_transition_anemia.png` — (c) Anemia 2×2 transition matrix

### Fig 9: Tuberculosis — ELA limitations for low-density single-basin diseases (Non-predictive)

**Caption:** (a) Temporal flow of basin composition among deceased TB patients.
Basin 1 retains ~68% of patients throughout the observation window, while
~32% remain in Non-Basin states at every time point — the highest stable
Non-Basin fraction across all 51 diseases. (b) PC sign pattern comparison
between the Basin state (S91) and the top 10 Non-Basin states. PC6 is the
sole axis distinguishing Basin from Non-Basin: PC6 = -1 in Basin, PC6 = +1
in all Non-Basin states, while PCs 1-5 and 7 vary. This suggests that
the landscape partitions TB patients along a single comorbidity axis (PC6)
rather than capturing disease severity. Cox analysis confirms HR > 1
(reversed direction), indicating that lower energy does not predict worse
outcomes in TB.

- `Fig09_tuberculosis.png` — Combined figure (assembled externally)
- `Fig09a_sankey_tuberculosis.png` — (a) TB Sankey (Basin 1 + Non-Basin)
- `Fig09b_pcpattern.png` — (b) PC sign pattern heatmap (PC6 highlighted)
- `Fig09c_pc6_loading.png` — (c) PC6 loading bar chart

## Supplementary Figures

- `FigS_atlas_hospitalization.png` — Atlas overview for hospitalization
- `FigS_jaccard_allblocks.png` — Jaccard clustering with 12 sub-cluster annotations (Fig 4 detail)
- `FigS_jaccard_basin_level.png` — Basin-level (93×93) Jaccard clustering; same-disease basins split apart
- `FigS_within_disease_jaccard.png` — Within-disease basin similarity dot plot with between-disease reference
- `FigS_code_traj_dementia_B{1,2}.png` — ICD-10 code frequency trajectory before death (Dementia, 2 basins)
- `FigS_code_traj_thyroid_B{1-4}.png` — ICD-10 code frequency trajectory before death (Thyroid, 4 basins)
- `FigS_code_traj_anemia_B{1,3}.png` — ICD-10 code frequency trajectory before death (Anemia, 2 basins)
- `FigS_code_traj_tuberculosis_B1.png` — ICD-10 code frequency trajectory before death (TB, 1 basin)
- `FigS_pattern_space.png` — Transition rate × asymmetry scatter (Table 1 backing)
- `FigS_death_hosp.png` — Death rate vs hospitalization rate across basins (Fig 7 generalization)
- `FigS_crossdisease_codes.png` — Top 10 cross-disease basin-discriminating codes (data-driven)
- `FigS_diabetes.png` — Diabetes (Convergent, overlaps Thyroid structure)
- `FigS_sepsis.png` — Sepsis (Structured non-predictive, HR >= 1)
- Per-disease plots for all 51 diseases in `plot/{disease}/`

## Supporting data files

- `Table1_pattern_classification.tsv` — 6-pattern classification for all 51 diseases

## Status

- [x] Fig 1: ChatGPT draft placed (Illustrator trace pending)
- [x] Fig 2: PCA combined figure placed
- [x] Fig 3: Atlas overview (per-basin majority, mortality)
- [x] Fig 4: Jaccard clustering (2-block dendrogram + heatmap)
- [x] Fig 5: Trajectory 3-panel (overview, basin start, non-basin start)
- [x] Fig 6: Thyroid bar chart + KM + disconnectivity + sankey (4 panels)
- [x] Fig 7: Dementia paired KM (mortality vs hospitalization)
- [x] Fig 8: Anemia transition trio (Symmetric → Convergent → Progressive)
- [ ] Fig 9: TB (Sankey + PC pattern — panels exist, combined figure TBD)
- [ ] Captions finalized for Fig 8-9
