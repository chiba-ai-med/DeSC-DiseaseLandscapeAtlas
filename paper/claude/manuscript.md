# Energy Landscape Analysis of 51 Disease Families Reveals Comorbidity Archetypes and Transition Dynamics from 11 Million Patient Claims

Koki Tsuyuzaki^1,2,\*, Katsuyuki Chida^3,4, Tetsuo Ishikawa^3,5,6,7, and Eiryo Kawakami^3,4,5,8,\*

^1 Data Science Core, Chiba University, 1-8-1 Inohana, Chuo, Chiba, Chiba 260-8670, JAPAN.
^2 Omics AI Research Team, Advanced General Intelligence in Science Program (AGIS), TRIP Headquarters, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan.
^3 Department of Artificial Intelligence Medicine, Graduate School of Medicine, Chiba University, 1-8-1 Inohana, Chuo, Chiba, Chiba 260-8670, JAPAN.
^4 Predictive Medicine Special Project (PMSP), RIKEN Center for Integrative Medical Sciences (IMS), RIKEN, Yokohama, Kanagawa, JAPAN.
^5 Division of Applied Mathematical Science, RIKEN Center for Interdisciplinary Theoretical and Mathematical Sciences (iTHEMS), RIKEN, Wako, Saitama, JAPAN.
^6 Department of Extended Intelligence for Medicine, The Ishii-Ishibashi Laboratory, Keio University School of Medicine, 35 Shinanomachi, Shinjuku-ku, Tokyo, 160-8582, Japan.
^7 Collective Intelligence Research Laboratory, Graduate School of Arts and Sciences, The University of Tokyo, 3-8-1 Komaba, Meguro-ku, Tokyo, 153-8902, Japan.
^8 Department of Advanced Artificial Intelligence Medicine, Graduate School of Medicine, The University of Osaka, 2-2 Yamada-oka, Suita, Osaka 565-0871, Japan.

\*Corresponding authors: Koki Tsuyuzaki (koki.tsuyuzaki@gmail.com) and Eiryo Kawakami (eiryo.kawakami@chiba-u.jp)

---

## Summary

Comorbidity patterns shape clinical trajectories, yet systematic cross-disease comparisons remain rare. We apply energy landscape analysis — an Ising model framework from statistical physics — to administrative claims from approximately 11 million individuals in the Japanese DeSC database, constructing comorbidity landscapes for 51 disease families. For each family, a binary comorbidity matrix is reduced to seven dimensions by principal component analysis, discretized into 128 states, and analyzed to identify 1–4 stable basins per disease. In 45 of 51 families, lower energy is associated with higher mortality risk in at least half of basins. Among deceased patients in single-basin diseases, first and last basin assignments matched in a median of 97% of cases, indicating that basins function as persistent attractors. Cross-disease comparison of basin comorbidity codes recovers clinically coherent disease clusters without supervision. We classify the 51 families into six behavioral patterns, providing a descriptive atlas of comorbidity organization across diseases.

**Keywords:** energy landscape analysis, comorbidity, claims data, disease phenotyping, unsupervised learning, Ising model

---

## 1 Introduction

Administrative claims databases capture diagnoses, procedures, and prescriptions for millions of individuals over extended periods, providing an increasingly important resource for studying disease at the population level (Schneeweiss and Avorn, 2005). While individual comorbidities and their pairwise associations have been extensively catalogued (Barnett et al., 2012), the higher-order structure of comorbidity — how multiple conditions co-occur, form stable configurations, and evolve over time — has received less attention. Network approaches have begun to map disease–disease relationships at population scale (Hidalgo et al., 2009; Jensen et al., 2014), but methods that identify discrete comorbidity states and their stability remain underexplored. Understanding this structure matters because comorbidity burden is a major driver of healthcare utilization, mortality, and treatment complexity (Salisbury et al., 2011).

Energy landscape analysis (ELA) offers a framework for discovering stable states in high-dimensional binary data. Originally developed in statistical physics to describe spin systems and applied to protein folding (Bryngelson et al., 1995), the approach constructs an energy function over discrete states such that frequently occupied configurations receive low energy. Local minima of this energy surface define "basins of attraction" — groups of states separated by energy barriers. The maximum-entropy (Ising model) formulation, which infers pairwise couplings from observed state frequencies (Schneidman et al., 2006), has been applied to neural population activity (Watanabe et al., 2014; Ezaki et al., 2017) and gene expression states (Sasai and Wolynes, 2003), among other domains (for a tutorial, see Masuda et al., 2025).

Here, we present the Disease Landscape Atlas, a systematic application of ELA to 51 disease families defined by ICD-10 code prefixes in the DeSC Japanese claims database. By constructing separate landscapes for each disease family — excluding the index disease codes from the feature set — the resulting basins reflect disease-specific comorbidity patterns rather than population-level prevalence. The cohorts collectively span approximately 11 million individuals observed between April 2014 and July 2021, with patients potentially contributing to multiple disease families. For each family, we construct a binary comorbidity matrix from 6-month observation windows, apply PCA to reduce dimensionality to 7 components, binarize the scores, and estimate the energy landscape using the Ising model.

Our analysis addresses three questions. First, does the direction of the energy–outcome association hold broadly across diseases, and for which disease families does it fail? Second, do the comorbidity codes defining basins reflect known clinical relationships, and does cross-disease comparison of basin content recover established disease groupings? Third, how stable are basin assignments over time, and do different diseases exhibit distinct patterns of inter-basin transitions?

---

## 2 Results

### 2.1 Analysis design and disease-family-specific landscape construction

We defined 51 disease families by ICD-10 code prefixes (Figure 1; Methods), spanning cardiovascular, neoplastic, neuropsychiatric, respiratory, metabolic, musculoskeletal, and infectious categories. For each family, an index cohort was constructed by identifying all individuals with at least one claim matching the family's defining ICD-10 codes in the DeSC database. Disease family cohorts ranged from 47,332 individuals (inflammatory bowel disease) to 4,514,492 individuals (diabetes mellitus) (Table S1). Because disease families are defined by non-exclusive ICD-10 prefixes and individuals can carry multiple diagnoses, patients may contribute to more than one family; the cohort sizes therefore cannot be summed to obtain the total number of unique individuals.

For each cohort, we constructed a binary comorbidity matrix in which rows represent patient–observation window pairs (6-month windows advanced in 1-month steps, yielding up to 83 windows per patient) and columns represent ICD-10 codes, excluding the index disease codes themselves. PCA reduced this matrix to 7 principal components per disease, and scores were binarized (positive → +1, non-positive → −1), yielding 2^7 = 128 discrete spin states per disease family. The Ising model energy function was estimated from the observed state frequencies using the Landscaper software, and local energy minima were identified as basins of attraction. Across the 51 families, the number of basins ranged from 1 to 4 (17 families with 1 basin, 26 with 2, 6 with 3, and 2 with 4).

### 2.2 PCA captures age and sex gradients across disease families

Although 7 principal components captured a modest fraction of total variance (median 12.5% across the 51 families as displayed in Figure 2c), they consistently captured biologically interpretable axes. PC1 correlated with patient age across nearly all disease families (median |r| = 0.41; Figure 2d), as illustrated by the diabetes cohort where older patients cluster at one pole of PC1 (Figure 2a). Sex separation was distributed across different PCs depending on the disease; in the urinary tract infection cohort, for example, PC7 separated female and male patients (|r| = 0.45; Figure 2b). The seven PCs were selected as a practical constraint: ELA requires enumerating 2^dim states, making 7 dimensions the practical maximum (128 states) for the Landscaper software. We do not claim that 7 components are optimal or capture sufficient variance for general-purpose dimensionality reduction; rather, they provide a tractable state space for landscape construction while preserving major clinical gradients.

### 2.3 Energy is associated with mortality and hospitalization across disease families

To assess whether the energy landscape captures clinically relevant structure, we tested the association between within-basin standardized energy and all-cause mortality using Cox proportional hazards models (adjusted for age and sex) within each basin separately (Figure 3; Methods). A disease family was classified as satisfying the per-basin majority criterion if HR < 1 (lower energy associated with higher hazard) held in at least half of its basins — including the case of exactly half. In 45 of 51 disease families, this criterion was satisfied for mortality. An analogous analysis using disease-specific hospitalization (primary diagnosis matching the index disease's ICD-10 codes; see Methods) as the endpoint yielded a consistent direction in 42 of 51 families (Figure S1).

The six disease families that did not satisfy the mortality criterion were hip fracture (1 basin, 0/1 HR < 1), liver cancer (3 basins, 0/3), lung cancer (1 basin, 0/1), pancreatic cancer (1 basin, 0/1), sepsis (3 basins, 0/3), and tuberculosis (1 basin, 0/1). These diseases tended to cluster in the region of low observation density, high terminality, or both (Figure 3), consistent with the hypothesis that ELA requires sufficient longitudinal comorbidity information to construct severity-aligned landscapes.

We emphasize that this criterion is a directional consistency check, not a predictive performance evaluation. The large sample sizes (up to 4.5 million patients per family) mean that even small and clinically negligible effect sizes can reach statistical significance; the hazard ratios should be interpreted as indicating the direction and approximate magnitude of the energy–outcome association, not as evidence of predictive utility.

### 2.4 Basin comorbidity codes recover clinically coherent disease clusters

To assess whether ELA basins capture clinically meaningful comorbidity structure, we compared the ICD-10 code compositions of basins across disease families. For each basin, we selected the 30 most frequent ICD-10 codes ranked by within-basin observation-window frequency (freq_in). For each pair of diseases, we computed the maximum Jaccard similarity between any two of their respective basins' code sets. This analysis excluded schizophrenia, yielding 50 disease families (Figure 4; Methods).

Hierarchical clustering of the resulting distance matrix (1 − maximum Jaccard; average linkage with optimal leaf ordering) revealed two major groups. An Acute/Frailty cluster contained 3 diseases: acute renal failure, aspiration pneumonia, and sepsis. The remaining 47 diseases formed a broad Chronic/Lifestyle cluster. Within the Chronic/Lifestyle cluster, clinically coherent sub-clusters emerged without supervision (Figure S2): neurodegenerative diseases (dementia, Parkinson disease, epilepsy), a frailty-adjacent group (delirium, pneumonia, malnutrition, hip fracture) positioned nearest the Acute/Frailty cluster, prostate diseases (benign prostatic hyperplasia, prostate cancer), pulmonary fibrosis diseases (interstitial lung disease, lung cancer), musculoskeletal diseases (osteoarthritis, osteoporosis), and a hepatic group (cirrhosis/fibrosis, liver cancer), among others.

These groupings, derived solely from the frequency of comorbidity codes within ELA basins, recapitulate known disease relationships. The separation of acute/frailty conditions from chronic/lifestyle conditions, and the sub-grouping of organ-system-related diseases, suggest that the basins capture clinically meaningful comorbidity structure. However, this is internal consistency within the same dataset, not external validation. At the individual basin level, same-disease basins did not always cluster together (Figure S3), indicating that different basins within a single disease can harbor comorbidity profiles more similar to basins of other diseases than to each other.

### 2.5 Longitudinal basin transitions and six behavioral patterns

To characterize the temporal stability of basin assignments, we tracked the basin membership of deceased patients across observation windows preceding death (Figure 5). Three observations emerged.

First, basins function as persistent attractors. Among deceased patients whose first observed state fell within a basin, a median of 97% had the same basin assignment at their first and last observation in single-basin diseases (17 families; Figure 5b). This first-to-last match rate permits intermediate departures and returns — it does not require continuous residence in the same basin. In diseases with 3–4 basins, the first-to-last match rate was approximately 69%, with 29–31% moving to a different basin by the final observation. Exit from basins to non-basin states was rare (0.4–3%). Second, non-basin states are transient: among deceased patients starting outside a basin, a median of 68–75% ended in a basin (Figure 5c). Third, the fraction of patients remaining entirely outside basins throughout follow-up (Stable Non-Basin in Figure 5a) was consistently below 2%. These statistics are medians across diseases within each basin-count group, computed from deceased patients only.

To organize the diversity of landscape behaviors across the 51 families, we defined a descriptive classification based on four axes: number of basins, direction of the energy–mortality association (per-basin majority criterion), inter-basin transition rate per step among deceased patients, and transition asymmetry (Table 1; Figure S4). This classification is descriptive rather than algorithmic: it was defined post hoc to summarize observed patterns, not discovered automatically or validated on independent data.

**Homogeneous** (13 families): Single basin, energy associated with mortality in the expected direction. These diseases — including COPD, stroke, ischemic heart disease, and colorectal cancer — have a single stable comorbidity configuration with a clear severity gradient within it.

**Reversed** (4 families): Single basin, reversed energy–mortality association. Hip fracture, lung cancer, pancreatic cancer, and tuberculosis have a single basin in which lower energy is not associated with worse outcomes.

**Convergent** (15 families): Multiple basins with expected energy–mortality direction, moderate-to-high transition rates (5–11% per step), and moderate transition asymmetry (0.3–0.6). These diseases — including thyroid disorder, diabetes, and epilepsy — exhibit net directional flow between basins.

**Symmetric** (11 families): Multiple basins with expected energy–mortality direction, moderate-to-high transition rates, but low asymmetry (<0.3). Dementia, atrial fibrillation, hypertension, CKD, and depression are examples: basins exchange patients at roughly equal rates.

**Progressive** (6 families): Multiple basins with expected energy–mortality direction, low transition rates (<5% per step), and high asymmetry (>0.6). Anemia, malignancy, and hematologic malignancy typify this pattern: inter-basin transitions are infrequent but strongly directional when they occur.

**Structured reversed** (2 families): Multiple basins, reversed energy–mortality association. Liver cancer and sepsis exhibit structured basin organization but the expected energy–outcome direction does not hold.

The thresholds for transition rate and asymmetry summarize observed variation in the data and are not derived from optimization criteria. Some diseases fall near the boundaries (e.g., peptic ulcer has a transition rate of 6.2% but is classified as Progressive; delirium has 4.6% but is classified as Convergent; heart failure has asymmetry 0.32, near the Symmetric/Convergent boundary).

### 2.6 Thyroid disorder: a Convergent landscape with clinically interpretable basins

The thyroid disorder cohort (n = 2,026,922; ICD-10 E00–E07, covering the full spectrum of thyroid diseases) yielded four basins, exemplifying the Convergent pattern (Figure 6). Basin-specific comorbidity profiles, ranked by the log-odds difference score (Methods), revealed clinically interpretable structure.

Basin 3 (n = 1,002,031; mean age 66.8) was characterized by codes related to chronic kidney disease and its complications: CKD stage 5 (N185; appearing in 1.3% of observation windows), hyperkalemia (E875; 3.2%), renal failure (N19; 6.1%), secondary hyperparathyroidism (E211; 1.7%), and chronic disease anemia (D638; 5.0%). Constipation (45.6% of observation windows) and xerosis cutis (21.1%) were also frequent in this basin (Figure 6a). The co-occurrence of thyroid disorder and renal codes is consistent with the known clinical association between thyroid dysfunction and CKD (Iglesias and Díez, 2009), although the cross-sectional design cannot establish directionality. Basin 3 had the lowest mortality among the four basins (39.4%).

Basin 4 (n = 135,495; mean age 76.0) was characterized by dental codes (chronic periodontitis, 81.8%; dental caries, 35.4%) along with hypertension (78.0%) and constipation (61.0%). Basin 2 (n = 108,187; mean age 74.7) had the highest mortality (49.4%), yielding a 10.0-percentage-point mortality gap relative to Basin 3 (Figure 6b). The disconnectivity graph showed that Basin 1 and Basin 3 have the lowest connection energy (E = −0.68), meaning these basins merge first in the hierarchical structure (Figure 6c). Temporal basin composition among deceased patients was largely stable across observation windows preceding death (Figure 6d).

The high frequency of dental codes in Basin 4 may reflect a population with regular healthcare contact; it does not establish that dental care is protective. Similarly, the co-occurrence patterns describe associations, not causal relationships.

### 2.7 Dementia: Symmetric pattern with mortality–hospitalization divergence

The dementia cohort (n = 933,285; ICD-10 F00–F03, G30) yielded two basins with a Symmetric transition pattern: Basin 1→2 transitions occurred at 5.7% per step and Basin 2→1 at 5.9% per step among deceased patients, indicating near-equal exchange (Figure 7; Figure 8a).

Despite this symmetry, the two basins showed different comorbidity profiles and hospitalization patterns. Basin 1 (n = 498,636; mean age 83.0) was characterized by constipation (79%), xerosis cutis (39%), dermatitis (38%), and heart failure (39%). Basin 2 (n = 434,649; mean age 80.8) was characterized by hypertension (78%), diabetes (48%), and hyperlipidemia (48%). The temporal evolution of these comorbidity profiles among deceased patients is shown in Figure S9.

Mortality was similar between the two basins (Basin 1: 49.3%; Basin 2: 53.9%; Figure 7a, left panel), with visually overlapping Kaplan–Meier curves despite a statistically significant log-rank test (p < 0.0001, driven by the large sample size). Hospitalization rates, however, diverged substantially: 77.6% in Basin 1 versus 57.2% in Basin 2 (Figure 7a, right panel). The hospitalization rate here refers to the proportion of patients with at least one admission of any cause during the observation period (all-cause hospitalization), which differs from the disease-specific hospitalization endpoint used in the per-basin majority criterion (§2.2). The Kaplan–Meier hospitalization panel in Figure 7a shows time to first admission. This mortality–hospitalization divergence suggests that the energy landscape captures dimensions of disease burden that are not fully aligned with terminal prognosis.

The Symmetric classification indicates that the two comorbidity profiles exchange patients at approximately equal rates. Whether these represent biologically distinct subtypes or different presentations of the same process cannot be determined from claims data alone.

### 2.8 Anemia: Progressive pattern with asymmetric transitions

Anemia exemplifies the Progressive pattern, with high transition asymmetry (Figure 8c). The anemia cohort (n = 2,251,643) had two basins: Basin 1 (n = 2,167,645) was large and characterized by CKD and dialysis-related codes (energy minimum = −0.990), while Basin 3 (n = 31,362) was small and characterized by dental and mild comorbidity codes (energy minimum = −0.100).

Among deceased patients, the conditional probability of transitioning from Basin 3 to Basin 1 at each step was 17.7%, whereas the probability of Basin 1 to Basin 3 was 0.5% — a ratio of 36:1 in conditional transition probabilities. However, because Basin 1 is much larger (69× more patients), the absolute numbers of transitioning patients were comparable: 71,040 (Basin 3→1) versus 70,562 (Basin 1→3) transition events. The asymmetry in conditional probabilities therefore reflects the large difference in basin sizes rather than a net unidirectional flow of patients. The overall transition rate was low (1.0% per step), and the mean transition asymmetry was 0.95 — the highest across all multi-basin diseases.

To illustrate how transition dynamics vary across ELA behavioral patterns, Figure 8 compares the step-by-step transition matrices of three representative diseases: dementia (Symmetric, Figure 8a), thyroid disorder (Convergent, Figure 8b), and anemia (Progressive, Figure 8c). In dementia, the two basins exchange patients at near-equal rates (5.7% vs. 5.9%). In thyroid disorder, moderate asymmetry is apparent: Basin 2→1 (8.9%) exceeds Basin 1→2 (1.9%). In anemia, the conditional probability asymmetry is extreme (36:1) but the absolute transition counts are balanced.

### 2.9 Tuberculosis: Reversed pattern and ELA limitations

Tuberculosis illustrates the limitations of ELA for diseases with low observation density and acute presentation. The TB cohort (n = 424,627) yielded a single basin (Basin 1) with a reversed energy–mortality association (HR = 1.069; 95% CI 1.063–1.074), meaning that lower energy was associated with lower, not higher, mortality (Figure 9).

Among deceased TB patients, 82.8% had trajectories classified as stable-in-basin (remaining in Basin 1 across all observed windows), while 3.3% remained outside Basin 1 throughout (trajectory_classification.tsv). The fraction of patients in non-basin states at any given time point is a separate quantity from the fraction who never enter a basin; Figure 9a shows the time-point composition, where the non-basin share visible at each time point should not be read as the percentage who remain in non-basin states permanently.

Analysis of the PC sign patterns distinguishing the basin state from the top-10 most-frequent non-basin states showed that PC6 was the sole consistently discriminating axis: the basin state had PC6 = −1, while all 10 examined non-basin states had PC6 = +1, with PCs 1–5 and 7 varying freely (Figure 9b). This comparison is limited to the selected states and does not establish that PC6 is the only discriminant across all 128 states.

The reversal of the energy–outcome direction in TB, combined with single-axis partitioning, exemplifies why certain diseases fail the per-basin majority criterion. TB patients may have relatively sparse claims records and rapid disease trajectories, providing insufficient longitudinal comorbidity information for landscape construction.

---

## 3 Discussion

We have presented the Disease Landscape Atlas, a systematic application of energy landscape analysis to 51 disease families in a Japanese claims database of approximately 11 million individuals. The atlas reveals that ELA basins capture clinically interpretable comorbidity structure in the majority of disease families and that this structure can be organized into six behavioral patterns.

**Basins as persistent comorbidity configurations.** Among deceased patients, basin assignments show high first-to-last concordance: a median of 97% in single-basin diseases and 69% in diseases with 3–4 basins (Figure 5b). Non-basin states are transient, with 68–75% of patients starting outside basins ending in one. This temporal persistence supports the physical analogy underlying ELA — that basins represent locally stable states in a comorbidity energy landscape. The degree of inter-basin exchange varies by disease, from the balanced exchange of dementia's Symmetric pattern to the highly asymmetric conditional transition probabilities of anemia's Progressive pattern. These observations are consistent with, but do not prove, specific clinical trajectories: in anemia, for example, the asymmetric transition probabilities are consistent with progressive renal decline driving chronic anemia through erythropoietin deficiency (Babitt and Lin, 2012), but the observational design cannot establish this mechanism.

**Cross-disease validity of basin structure.** The Jaccard clustering analysis (Figure 4) shows that basin comorbidity profiles reflect recognizable clinical relationships. Clinically coherent clusters — neurodegenerative diseases, prostate diseases, musculoskeletal diseases, and a tight acute/frailty core (acute renal failure, aspiration pneumonia, sepsis) — emerged without supervision, echoing disease groupings identified by network approaches using different data and methods (Hidalgo et al., 2009; Jensen et al., 2014). Frailty-associated conditions (delirium, pneumonia, malnutrition, hip fracture) fell within the Chronic/Lifestyle cluster but at its boundary nearest the Acute/Frailty group, consistent with their dual chronic-and-acute character. However, this internal consistency does not constitute external validation, as the same data that defined the basins also generated the cluster structure.

**Energy–outcome association and its limitations.** The per-basin majority criterion was satisfied in 45 of 51 families for mortality and 42 of 51 for hospitalization, suggesting that the energy–outcome direction is broadly consistent. The six families failing the mortality criterion include single-basin diseases with acute presentation or high lethality (hip fracture, lung cancer, pancreatic cancer, tuberculosis) and multi-basin diseases where no basin showed the expected direction (liver cancer, sepsis). This pattern is consistent with the interpretation that ELA requires sufficient longitudinal comorbidity information to construct severity-aligned landscapes.

**Mortality–hospitalization divergence.** The dementia case study (Figure 7) highlights that basins with similar mortality can differ substantially in all-cause hospitalization rates: 77.6% versus 57.2% despite a mortality difference of only 4.6 percentage points. Whether this reflects genuine differences in disease burden, care-seeking behavior, or healthcare system factors cannot be determined from claims data alone. A systematic comparison across multi-basin diseases (Figure S5) showed that mortality–hospitalization alignment varies widely.

**The six-pattern classification.** The categories provide a descriptive vocabulary for the behavioral diversity of comorbidity landscapes. We emphasize that this classification was defined post hoc and should not be treated as a validated taxonomy. The boundaries between categories (e.g., asymmetry thresholds) were chosen to summarize observed variation, not derived from external criteria. Some diseases fall near category boundaries. Future work could assess whether these categories are robust to different analytical parameters or data sources.

**Limitations.** Several limitations warrant consideration. First, ICD-10 codes from administrative claims are recorded for billing purposes and may not reflect clinical reality; coding practices vary across institutions and over time (O'Malley et al., 2005). Second, the overlapping observation windows (6-month windows sliding by 1 month) share 5 of 6 months of information, which may inflate apparent temporal stability. Third, the 7-dimensional PCA reduction captures a median of 12.5% of total variance (13.3% by eigenvalue calculation; Figure 2c shows the patient-representative estimate). The 128 discrete states are a coarse approximation of the comorbidity space. Fourth, disease families are defined by overlapping ICD-10 prefixes, so patients can appear in multiple families. Fifth, the patient-level representative record for Cox analyses selects the observation window with the maximum row identifier (max old_row_id), whose correspondence to calendar ordering has not been independently verified. This record may incorporate comorbidity information from the patient's observation history and should not be interpreted as a baseline measurement. Sixth, we did not perform proportional hazards assumption testing, competing risk analysis, or multiple testing correction across the 51 disease families. We report nominal p-values and focus on effect direction and magnitude rather than statistical significance. Seventh, the longitudinal trajectory and transition analyses are restricted to deceased patients; whether the same patterns hold in the full population is not established. Eighth, the energy landscape is estimated from the static cross-sectional distribution of states and does not directly model temporal dynamics. The connection energies in the disconnectivity graph describe the energy surface, not measured transition rates.

**Outlook.** The Disease Landscape Atlas demonstrates that ELA can extract stable, interpretable, and cross-disease-comparable comorbidity structures from claims data at scale. Potential applications include identifying patient subgroups for stratified analysis, generating hypotheses about comorbidity-driven disease trajectories, and comparing comorbidity architecture across populations or healthcare systems. Validation in independent cohorts and formal assessment of the relationship between energy barriers and transition dynamics represent important next steps.

---

## 4 Methods

### 4.1 Data source and construction of disease-specific matrices

The analysis used administrative claims data from the DeSC database, a Japanese claims and health checkup database, covering the period from April 2014 to July 2021. The database contains diagnosis codes (ICD-10, 2013 Japanese modification), demographic information (birth year-month, sex, insurer type), observation period (start and end of coverage), death flag, and hospitalization records.

We defined 51 disease families based on ICD-10 code prefixes (Supplementary Table S1). Each family was specified by one or more ICD-10 prefixes; for example, the dementia family included F00–F03 and G30, while the diabetes family included E10–E14. For each disease family, an index cohort was defined as all individuals with at least one claims record matching the family's ICD-10 prefix codes during the observation period. Observation windows of 6 months' duration were constructed with a sliding step of 1 month, yielding up to 83 windows per patient across the full observation period.

A binary comorbidity matrix was constructed for each cohort, with rows representing patient–observation window pairs and columns representing ICD-10 diagnosis codes. The matrix entry was 1 if the code appeared in at least one claim within that observation window and 0 otherwise. The ICD-10 codes defining the index disease were excluded from the column set for each family, so that the resulting landscape reflects comorbidity patterns rather than the index disease itself. After construction, the matrix was stored in coordinate (COO) format with zstd compression.

### 4.2 Principal component analysis

For each disease family, the comorbidity matrix was subjected to exact, out-of-core PCA using OnlinePCA.jl (Tsuyuzaki, 2026; Tsuyuzaki et al., 2020) with a chunk size of 300,000 rows. The first 7 principal components were retained for all families. The choice of 7 dimensions was driven by the Landscaper software's requirement to enumerate all 2^dim states; 7 dimensions yield 128 states. The median cumulative variance explained by 7 PCs was 12.5% as displayed in Figure 2c (computed from patient-representative scores; eigenvalue-based calculation yields 13.3%, the difference reflecting deduplication to one record per patient).

### 4.3 Binarization of principal component scores

Each PCA score was binarized using a zero threshold: values greater than 0 were assigned B = 1, and values less than or equal to 0 were assigned B = 0. The corresponding Ising spin variable was defined as S = 2B − 1, yielding S ∈ {−1, +1} for each of the 7 components. Each patient–observation window was thereby assigned one of 2^7 = 128 discrete states. The Landscaper software expects a sparse binary matrix where non-zero entries indicate S = +1 and absent entries indicate S = −1.

### 4.4 Energy landscape analysis using Landscaper

Energy landscape estimation was performed using Landscaper (https://github.com/chiba-ai-med/Landscaper), which fits an Ising model to the observed state frequencies and computes the energy E(σ) for each of the 128 states σ. States with lower energy correspond to higher observation probability. Local energy minima were identified, and each state was assigned to the basin of the nearest local minimum under steepest-descent dynamics. The SubGraph.tsv output records basin labels for states that belong to a non-trivial basin. In 47 of 51 diseases, SubGraph.tsv excludes isolated local minima (states with no lower-energy neighbor in the StatusNetwork), yielding fewer than 128 rows; in the remaining 4 diseases (dementia, delirium, diabetes, sepsis), all 128 states are included. Basin.tsv records the state index of each local minimum. The disconnectivity graph shows the connection energies at which basin pairs merge. Landscaper was run with a single thread per disease.

### 4.5 Basin assignment and clinical outcome analysis

**Patient-level representative records.** Because each patient contributed multiple observation windows, a single representative record per patient was required for clinical outcome analysis. For each patient, the observation window with the maximum row identifier (max old_row_id) was selected. This record's state, basin assignment, and energy value served as the patient-level summary. The relationship between row identifier ordering and calendar time has not been independently verified; consequently, this record should not be assumed to correspond to the latest or earliest observation window without further confirmation.

**Cox proportional hazards models.** For each disease family and each basin, a Cox proportional hazards model was fitted:

> h(t) = h₀(t) × exp(β₁ E_z + β₂ age + β₃ sex)

where E_z is the energy value standardized (Z-score) within each basin, age is calculated from birth year-month to coverage start date, and sex is coded as a factor. Two endpoints were analyzed separately:

*Mortality.* The time variable is defined as the duration of the patient's insurance coverage period (observable_end_ym − observable_start_ym, converted to years). The event indicator is shibou_flg_riyouka, a derived death flag in the DeSC claims database. This flag is set for a broader population than the original death flag (shibou_flg) in the source data; the derivation criteria are defined by the database provider. The code does not compute a separate death date or take a minimum of death date and coverage end date — it uses the coverage end month (observable_end_ym) as the time endpoint for all patients, with shibou_flg_riyouka as the event indicator. For deceased patients, the relationship between observable_end_ym and the actual time of death is determined by the database construction and has not been independently verified.

*Disease-specific hospitalization.* The time origin is likewise the start of coverage. The event is first admission where the DPC primary diagnosis matches the index disease family's ICD-10 codes. For patients with a qualifying admission, the time variable is the interval from coverage start to admission. For patients without a qualifying admission, the time variable is the full coverage duration (observable_end_ym − observable_start_ym); death is not explicitly modeled as a competing risk or separate censoring event in this analysis.

The case-study hospitalization rates reported for individual diseases (e.g., dementia Basin 1: 77.6%, Basin 2: 57.2% in Figure 7) use all-cause hospitalization (any admission, regardless of diagnosis) from a separate population-wide hospitalization register, which differs from the disease-specific endpoint used in the per-basin majority criterion.

**Per-basin majority criterion.** A disease family was classified as satisfying the criterion if HR < 1 held in at least half of its basins. This includes ties: for a disease with 2 basins, 1 of 2 suffices; for 4 basins, 2 of 4 suffices.

**Kaplan–Meier curves.** Kaplan–Meier survival curves were constructed for inter-basin comparison, stratified by basin assignment from the representative record. For visualization, energy tertile groups (low, medium, high) within each basin were used to illustrate intra-basin variation.

### 4.6 Comorbidity profiles and cross-disease similarity

**Log-odds difference scoring.** To characterize each basin's comorbidity profile, we computed a log-odds difference score for each ICD-10 code j in basin b:

> score_{bj} = logit((c_{bj} + α) / (M_b + 2α)) − logit((c_{all,j} + α) / (M_{all} + 2α))

where c_{bj} is the number of observation windows in basin b containing code j, M_b is the total count of all code occurrences across observation windows in basin b, c_{all,j} and M_{all} are the corresponding totals across all basins (including basin b itself), α = 1 is a smoothing parameter, and logit(p) = ln(p / (1 − p)). Codes were filtered to those with within-basin observation-window frequency > 1% (freq_in > 0.01) and at least 100 occurrences (n_in ≥ 100). The denominator for freq_in is the number of observation windows assigned to the basin, not the number of unique patients.

**Jaccard similarity and cross-disease clustering.** For each basin, the top 30 ICD-10 codes were selected by within-basin observation-window frequency (freq_in), not by log-odds difference score, and without additional filters beyond the top-30 cutoff. For each pair of disease families, the maximum Jaccard index between any two of their basins' top-30 code sets was computed. The resulting similarity matrix was converted to a distance matrix (1 − Jaccard) and subjected to average-linkage hierarchical clustering with optimal leaf ordering. Schizophrenia was excluded, leaving 50 disease families.

### 4.7 Longitudinal trajectory and transition analysis

**Trajectory-to-death analysis.** For each deceased patient in each disease family, we tracked the patient's state, energy, and basin assignment across retained observation windows, indexed by ordinal position counting backward from the last retained observation (periods_before_death = 0 at the final retained window). The ordinal positions reflect the sequence of available records, not calendar-time intervals; gaps in coverage or excluded windows are not explicitly represented. Up to 20 positions before death were retained per patient.

**Step-by-step transition matrices.** Using the deceased-patient trajectory data, we computed transition rates between basin assignments at consecutive retained positions. For each pair of consecutive positions (t and t+1 in ordinal sequence), the patient's basin label at each position determined a transition event. Because observation windows are constructed with a 1-month sliding step, we report transition rates "per step" throughout, referring to the interval between consecutive retained observation positions. The calendar-time gap between consecutive retained positions has not been independently verified and may vary across patients; "per step" should not be read as "per calendar month." The transition rate from basin i to basin j was defined as the count of i→j transitions divided by the total transitions departing basin i. Transition asymmetry for each basin pair (i, j) was defined as |rate(i→j) − rate(j→i)| / (rate(i→j) + rate(j→i)). For diseases with more than 2 basins, the mean asymmetry across all basin pairs was used. These analyses are restricted to deceased patients; results may differ in the full population.

**Figure 5 trajectory summaries.** Figure 5a shows the distribution of trajectory categories (Stable Basin, Stable Non-Basin, Moved) as medians across diseases within each basin-count group. Figures 5b and 5c compare the first and last observed basin assignments (endpoint comparison), which permits intermediate changes. These summaries are medians across diseases, not pooled patient-level statistics.

**ICD-10 code frequency trajectory heatmaps.** To visualize how comorbidity composition changes over time before death within each basin, we computed the observation-window-level frequency of each ICD-10 code at each ordinal position before death. For each basin in each disease family, the patient set consisted of all deceased patients assigned to that basin at each retained position (basin membership was re-evaluated at each time point, so a patient who changed basins contributed to different basins at different positions). The denominator for frequency was the number of unique patients with a retained observation at that position within the basin. The numerator was the number of observation windows in that basin–position combination containing the code (one count per observation window, regardless of how many times the code appeared). The top 30 ICD-10 codes were selected by total occurrence count across all positions within the basin. These heatmaps show changes in population-level code frequencies over time; they do not track individual patients' disease onset or measure pairwise comorbidity co-occurrence.

### 4.8 Software and reproducibility

The analysis pipeline was implemented using Snakemake (Mölder et al., 2021) with Singularity containers for reproducibility. Key software components were:

- **OnlinePCA.jl** (container: `ghcr.io/rikenbit/onlinepcajl:df65589`): exact out-of-core PCA (Tsuyuzaki, 2026; Tsuyuzaki et al., 2020)
- **Landscaper** (container: `ghcr.io/chiba-ai-med/landscaper:main`): Ising model energy landscape estimation
- **R** (survival, data.table, ggplot2 packages): Cox analysis, Kaplan–Meier estimation, visualization
- **Julia** (MatrixMarket.jl, SparseArrays.jl, CodecZstd.jl): sparse matrix I/O, PCA score processing, binarization

All workflow definitions are available in the project repositories: https://github.com/chiba-ai-med/DeSC-DiseaseLandscapeAtlas (atlas analysis) and https://github.com/chiba-ai-med/DeSC-investigation (data preprocessing). Configuration is specified in `config/config.yaml` and `config/disease_families.tsv`.

---

## 5 Resource availability

### Lead contact

Further information and requests for resources should be directed to Koki Tsuyuzaki (koki.tsuyuzaki@gmail.com).

### Materials availability

This study did not generate new reagents or materials.

### Data and code availability

- The analysis code and workflow definitions are available at https://github.com/chiba-ai-med/DeSC-DiseaseLandscapeAtlas (atlas analysis), https://github.com/chiba-ai-med/DeSC-investigation (data preprocessing), https://github.com/rikenbit/OnlinePCA.jl (PCA), and https://github.com/chiba-ai-med/Landscaper (energy landscape estimation).
- The DeSC database, which includes claims data collected during routine medical practice, is not publicly available. Access to the dataset requires a Data Sharing Agreement with DeSC Healthcare, Inc. However, data required to reproduce the analysis may be shared upon reasonable request and with the permission of DeSC Healthcare, Inc.
- Summary-level data supporting the findings of this study (atlas summary tables, basin profiles, transition matrices, and Jaccard similarity matrices) are available in the project repository.
- Any additional information required to reanalyze the data reported in this paper is available from the lead contact upon request.

---

## Acknowledgments

This work was supported by JSPS KAKENHI (grant numbers 23K11312, 24K02011, 24K21269, and 23K23933), JST A-STEP (grant number JPMJTR25UJ), and JST PRESTO (grant number JPMJPR1945) to K.T.

This study was approved as exempt from ethical review by the institutional review board, as it used only anonymized claims data with no individually identifiable information.

---

## Author contributions

K.T.: Conceptualization, Methodology, Software, Formal analysis, Visualization, Writing – original draft. K.C.: Data curation, Resources. T.I.: Methodology, Supervision. E.K.: Methodology, Supervision. All authors: Investigation, Writing – review & editing.

---

## Declaration of interests

This study was conducted under a joint research agreement with DeSC Healthcare, Inc., which provided access to the DeSC database. DeSC Healthcare, Inc. had no role in the study design, analysis, interpretation of the data, or the decision to submit the manuscript for publication. The authors have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

---

## Declaration of generative AI use

During the preparation of this work, the authors used generative AI tools as follows: Claude Code (Anthropic) for data analysis pipeline development and manuscript drafting, and ChatGPT (OpenAI) and Claude Code for manuscript editing and revision. The authors reviewed and edited all AI-generated content, take full responsibility for the content of this publication, and confirm that the final manuscript reflects their own scientific judgment.

---

## References

1. Babitt JL, Lin HY. Mechanisms of anemia in CKD. *Journal of the American Society of Nephrology*. 2012;23(10):1631–1634. doi:10.1681/ASN.2011111078. PMID: 22935483.

2. Barnett K, Mercer SW, Norbury M, Watt G, Wyke S, Guthrie B. Epidemiology of multimorbidity and implications for health care, research, and medical education: a cross-sectional study. *The Lancet*. 2012;380(9836):37–43. doi:10.1016/S0140-6736(12)60240-2. PMID: 22579043.

3. Bryngelson JD, Onuchic JN, Socci ND, Wolynes PG. Funnels, pathways, and the energy landscape of protein folding: a synthesis. *Proteins*. 1995;21(3):167–195. doi:10.1002/prot.340210302. PMID: 7784423.

4. Ezaki T, Watanabe T, Ohzeki M, Masuda N. Energy landscape analysis of neuroimaging data. *Philosophical Transactions of the Royal Society A*. 2017;375(2096):20160287. doi:10.1098/rsta.2016.0287.

5. Hidalgo CA, Blumm N, Barabási AL, Christakis NA. A dynamic network approach for the study of human phenotypes. *PLoS Computational Biology*. 2009;5(4):e1000353. doi:10.1371/journal.pcbi.1000353. PMID: 19360091.

6. Iglesias P, Díez JJ. Thyroid dysfunction and kidney disease. *European Journal of Endocrinology*. 2009;160(4):503–515. doi:10.1530/EJE-08-0837. PMID: 19095779.

7. Jensen AB, Moseley PL, Oprea TI, et al. Temporal disease trajectories condensed from population-wide registry data covering 6.2 million patients. *Nature Communications*. 2014;5:4022. doi:10.1038/ncomms5022. PMID: 24959948.

8. Masuda N, Islam S, Aung ST, Watanabe T. Energy landscape analysis based on the Ising model: tutorial review. *PLOS Complex Systems*. 2025;2(5):e0000039. doi:10.1371/journal.pcsy.0000039.

9. Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33. doi:10.12688/f1000research.29032.2.

10. O'Malley KJ, Cook KF, Price MD, Wildes KR, Hurdle JF, Ashton CM. Measuring diagnoses: ICD code accuracy. *Health Services Research*. 2005;40(5 Pt 2):1620–1639. doi:10.1111/j.1475-6773.2005.00444.x. PMID: 16178999.

11. Salisbury C, Johnson L, Purdy S, Valderas JM, Montgomery AA. Epidemiology and impact of multimorbidity in primary care: a retrospective cohort study. *British Journal of General Practice*. 2011;61(582):e12–e21. doi:10.3399/bjgp11X548929.

12. Sasai M, Wolynes PG. Stochastic gene expression as a many-body problem. *Proceedings of the National Academy of Sciences*. 2003;100(5):2374–2379. doi:10.1073/pnas.2627987100.

13. Schneeweiss S, Avorn J. A review of uses of health care utilization databases for epidemiologic research on therapeutics. *Journal of Clinical Epidemiology*. 2005;58(4):323–337. doi:10.1016/j.jclinepi.2004.10.012. PMID: 15862718.

14. Schneidman E, Berry MJ II, Segev R, Bialek W. Weak pairwise correlations imply strongly correlated network states in a neural population. *Nature*. 2006;440(7087):1007–1012. doi:10.1038/nature04701. PMID: 16625187.

15. Tsuyuzaki K. OnlinePCA.jl: A Julia Package for Out-of-core and Sparse Principal Component Analysis. *Journal of Open Source Software*. 2026;11(117):9343. doi:10.21105/joss.09343.

16. Tsuyuzaki K, Sato H, Sato K, Nikaido I. Benchmarking principal component analysis for large-scale single-cell RNA-sequencing. *Genome Biology*. 2020;21:9. doi:10.1186/s13059-019-1900-3. PMID: 31955711.

17. Watanabe T, Masuda N, Megumi F, Kanai R, Rees G. Energy landscape and dynamics of brain activity during human bistable perception. *Nature Communications*. 2014;5:4765. doi:10.1038/ncomms5765. PMID: 25163855.

---

## Figure legends

### Figure 1. Overview of the Disease Landscape Atlas analysis workflow.

(a) From the DeSC claims database, 51 disease-family cohorts are defined by ICD-10 index codes. For each cohort, a binary comorbidity matrix is constructed in which rows represent patient–observation window pairs and columns represent ICD-10 codes excluding the index disease. PCA reduces the matrix to 7 dimensions, scores are binarized to ±1, yielding 2^7 = 128 discrete states. Energy landscape analysis assigns each state an energy value and groups states into basins of attraction (1–4 per disease). (b) Within each basin, Cox proportional hazards models test whether lower energy is associated with higher mortality or hospitalization risk (HR < 1 per 1 SD increase in energy, adjusted for age and sex). The per-basin majority criterion (HR < 1 in at least half of basins) is satisfied in 45 of 51 families for mortality and 42 of 51 for hospitalization.

### Figure 2. PCA captures age and sex axes across disease families.

(a) Diabetes cohort projected onto PC1 vs. PC2, colored by age. PC1 correlates with age across most diseases (median |r| = 0.41). (b) Urinary tract infection cohort projected onto PC1 vs. PC7, colored by sex. PC7 separates female and male patients (|r| = 0.45). (c) Cumulative variance explained by the first 7 PCs for all 51 disease families (grey) and their median (red). Seven PCs capture a median of 12.5% of total variance. (d) Pearson correlation of each PC with age (left) and sex (right) across 51 disease families, sorted by PC1–age correlation. PC1 universally captures an age gradient; sex associations are distributed across different PCs per disease.

### Figure 3. Atlas of 51 disease families by observation density, death rate, and energy–mortality association.

Each point represents one disease family. The x-axis shows the mean number of claims records per patient (observation density), and the y-axis shows the cohort death rate during follow-up. Point color indicates whether energy is associated with mortality in the expected direction (HR < 1) in at least half of basins (blue) or not (red). The per-basin majority criterion is satisfied in 45 of 51 disease families.

### Figure 4. Disease clustering by basin code similarity recovers clinically valid disease groups.

Hierarchical clustering of 50 disease families (schizophrenia excluded) by basin code similarity. For each pair of diseases, the maximum Jaccard similarity between any two basins — using each basin's top-30 most frequent ICD-10 codes by within-basin observation-window frequency — is computed. Average-linkage clustering with optimal leaf ordering reveals two major groups: an Acute/Frailty cluster (3 diseases: acute renal failure, aspiration pneumonia, sepsis) and a broad Chronic/Lifestyle cluster (47 diseases). Within the Chronic/Lifestyle cluster, clinically coherent sub-clusters emerge without supervision, including a frailty-adjacent group at the boundary nearest the Acute/Frailty cluster.

### Figure 5. Basin assignments persist among deceased patients, with basin-count-dependent transition rates.

(a) Distribution of deceased patients' trajectory categories by number of basins. Each bar shows the median across diseases in that group. In 1-basin diseases (17 families), a median of 88% remain in the same basin throughout follow-up; in 2+ basin diseases, roughly half undergo state transitions. Stable Non-Basin accounts for less than 2%. (b) Among deceased patients starting in a basin, a median of 97% had the same basin at first and last observation in 1-basin diseases; this drops to approximately 69% in 3–4 basin diseases, where 29–31% move to a different basin. This is an endpoint comparison that permits intermediate departures and returns. (c) Among deceased patients starting outside a basin, a median of 68–75% ended in a basin. n = 24 for the 2-basin group, as 2 of 26 two-basin families had no non-basin-starting deceased patients.

### Figure 6. Thyroid disorder energy landscape with clinically interpretable basins (Convergent pattern).

(a) Top 10 comorbidity codes in Basin 3 (n = 1,002,031) and Basin 4 (n = 135,495), ranked by log-odds difference score. Bar color indicates within-basin observation-window frequency. Basin 3 is characterized by CKD-related codes including CKD stage 5, hyperkalemia, renal failure, secondary hyperparathyroidism, and chronic disease anemia. Basin 4 is characterized by dental codes (chronic periodontitis 81.8%, dental caries 35.4%). (b) Kaplan–Meier survival curves for all four basins. Basin 2 (death rate 49.4%, mean age 74.7) has the highest mortality, while Basin 3 (death rate 39.4%, mean age 66.8) has the lowest, a 10.0-percentage-point gap. (c) Disconnectivity graph. Basin 1 and Basin 3 merge at the lowest connection energy (E = −0.68), indicating energetic adjacency. (d) Temporal basin composition among deceased patients at successive observation positions before the last retained record.

### Figure 7. Dementia — similar mortality but divergent hospitalization between two basins (Symmetric pattern).

(a) Kaplan–Meier curves for mortality (left) and all-cause hospitalization (right) in the dementia cohort (n = 933,285), stratified by basin. The two basins show similar mortality (Basin 1: 49.3%, Basin 2: 53.9%) but different all-cause hospitalization patterns (Basin 1: 77.6%, Basin 2: 57.2%). Note that these all-cause admission rates differ from the disease-specific hospitalization endpoint used in the per-basin majority criterion (§2.2). Basin 1 (mean age 83.0) is associated with constipation, xerosis cutis, dermatitis, and heart failure codes. Basin 2 (mean age 80.8) is associated with hypertension, diabetes, and hyperlipidemia codes. Transition rates between basins among deceased patients are near-symmetric (5.7% vs. 5.9% per step). The hospitalization panel shows Kaplan–Meier curves for time to first all-cause admission, with censoring at the end of the observation period for patients without an admission.

### Figure 8. Transition dynamics across three ELA behavioral patterns.

Step-by-step basin transition matrices among deceased patients for three representative diseases, ordered by increasing asymmetry. (a) Dementia (Symmetric): Basin 1→2 = 5.7%, Basin 2→1 = 5.9%, near-symmetric exchange. (b) Thyroid disorder (Convergent): 4×4 matrix with moderate asymmetry; Basin 2→1 (8.9%) exceeds Basin 1→2 (1.9%). (c) Anemia (Progressive): Basin 3→1 = 17.7% vs. Basin 1→3 = 0.5%, a 36:1 ratio in conditional transition probabilities. The absolute transition counts are comparable (71,040 vs. 70,562) because Basin 1 is much larger than Basin 3.

### Figure 9. Tuberculosis — ELA limitations for low-density single-basin diseases (Reversed pattern).

(a) Temporal basin composition among deceased TB patients at successive observation positions. (b) PC sign pattern comparison between the basin state and the top-10 most-frequent non-basin states. PC6 is the sole axis that consistently differs: PC6 = −1 in the basin state vs. PC6 = +1 in all 10 examined non-basin states. This comparison is limited to the selected states. (c) PC6 loading bar chart showing the ICD-10 codes contributing to the discriminating axis. Cox analysis yields HR = 1.069 (95% CI: 1.063–1.074), indicating reversed energy–mortality direction.

---

```{=openxml}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="11906" w:h="16838" w:orient="portrait"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

## Table 1. Six-pattern classification of 51 disease families

| Pattern | Basins | Cox | Rate | Asym. | n | Diseases |
|---------|--------|-----|------|-------|---|----------|
| Homogeneous | 1 | Exp | — | — | 13 | Autoimmune conn. tissue, Breast ca., Cirrhosis, Colorectal ca., COPD, Gout, ILD, IHD, Liver disease, Metastatic ca., Obesity, Stomach ca., Stroke |
| Reversed | 1 | Rev | — | — | 4 | Hip fracture, Lung ca., Pancreatic ca., TB |
| Progressive | Multi | Exp | <5% | >0.6 | 6 | Anemia, Hem. malignancy, IBD, Malignancy, Peptic ulcer, UTI |
| Convergent | Multi | Exp | 5–11% | 0.3–0.6 | 15 | Asp. pneumonia, BPH, Delirium, Diabetes, Epilepsy, Heart failure, Malnutrition, OA, Osteoporosis, Pneumonia, RA, Schizophrenia, Sleep disorders, Thyroid, Viral hepatitis |
| Symmetric | Multi | Exp | 5–10% | <0.3 | 11 | AF, AKI, Anxiety, Asthma, CKD, Dementia, Depression, Dyslipidemia, HTN, Parkinson, Prostate ca. |
| Str. reversed | Multi | Rev | Any | Any | 2 | Liver ca., Sepsis |

Abbreviations: Exp = Expected direction (HR < 1 in ≥ half of basins; directional consistency, not predictive performance); Rev = Reversed; Rate = inter-basin transition rate per step; Asym. = mean pairwise transition asymmetry; ca. = cancer; conn. = connective; ILD = interstitial lung disease; IHD = ischemic heart disease; Hem. = hematologic; IBD = inflammatory bowel disease; UTI = urinary tract infection; Asp. = aspiration; AF = atrial fibrillation; BPH = benign prostatic hyperplasia; OA = osteoarthritis; RA = rheumatoid arthritis; AKI = acute renal failure; HTN = hypertension; Str. reversed = Structured reversed; TB = tuberculosis. Thresholds summarize observed variation; some diseases fall near the stated boundaries. Single-basin diseases have no transition metrics. This classification is descriptive and was not pre-registered or validated on independent data.

```{=openxml}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

---

## Supplementary figure legends

### Figure S1. Atlas overview for hospitalization.

Same format as Figure 3, but using disease-specific hospitalization (DPC primary diagnosis matching the index disease's ICD-10 codes) as the endpoint. The per-basin majority criterion (HR < 1 in ≥ half of basins) is satisfied in 42 of 51 disease families.

### Figure S2. Jaccard clustering with sub-cluster annotations.

Same clustering as Figure 4, with additional annotations highlighting clinically coherent sub-clusters within the Chronic/Lifestyle group: a frailty-adjacent group (delirium, pneumonia, malnutrition, hip fracture), neurodegenerative (dementia, Parkinson disease, epilepsy), prostate (benign prostatic hyperplasia, prostate cancer), pulmonary fibrosis (interstitial lung disease, lung cancer), musculoskeletal (osteoporosis, osteoarthritis), hepatic (cirrhosis/fibrosis, liver cancer), and others.

### Figure S3. Basin-level Jaccard clustering.

Hierarchical clustering of individual basins (across 50 diseases) by top-30 code Jaccard similarity. Basins from the same disease do not necessarily cluster together, indicating that within-disease basins can harbor comorbidity profiles more similar to basins of other diseases than to each other.

### Figure S4. Transition rate vs. asymmetry scatter (Table 1 backing).

Inter-basin transition rate per step (x-axis) vs. mean pairwise asymmetry (y-axis) for multi-basin disease families. Points are colored by Table 1 pattern classification.

### Figure S5. Death rate vs. hospitalization rate across basins.

Scatter plot of per-basin death rate (x-axis) vs. per-basin all-cause hospitalization rate (y-axis) across all basins in multi-basin disease families. Hospitalization rates here are all-cause (any admission).

### Figure S6. Top cross-disease basin-discriminating codes.

Bar chart showing the ICD-10 codes most frequently appearing in basins' top-30 lists across multiple disease families.

### Figure S7. Diabetes — Convergent pattern (supplementary case study).

Four-panel figure for the diabetes cohort (n = 4,514,492): (a) basin comorbidity profiles, (b) inter-basin Kaplan–Meier curves, (c) disconnectivity graph, (d) temporal basin flow.

### Figure S8. Sepsis — Structured reversed pattern (supplementary case study).

Sepsis cohort (n = 434,910) with 3 basins (Basin 1: 154,734; Basin 2: 233,467; Basin 3: 46,709). Despite structured basin organization, the energy–mortality direction is reversed (HR > 1 in all basins). Panels show: KM survival curves by basin, tag clouds of top comorbidity codes per basin, and trajectory heatmap before death.

### Figure S9. ICD-10 code frequency trajectories before death and before hospitalization, per basin.

Heatmaps showing the frequency of the top 30 ICD-10 codes (by total count within each basin) across ordinal positions before death (panels a) and before disease-specific hospitalization (panels b), for representative disease families (dementia, thyroid disorder, anemia, tuberculosis). Each cell shows the fraction of observation windows at that position within the basin that contained the code (denominator: number of unique patients at that position in the basin). The x-axis represents ordinal positions before the event (0 = last retained observation before death or hospitalization), not calendar time; the y-axis lists codes ordered by overall frequency. Basin membership is re-evaluated at each position, so patients who changed basins contribute to different basins at different positions. The hospitalization endpoint is disease-specific: admission where the DPC primary diagnosis matches the index disease family's ICD-10 codes. These heatmaps reflect population-level shifts in code prevalence and do not imply individual-level disease onset ordering or pairwise comorbidity co-occurrence.
