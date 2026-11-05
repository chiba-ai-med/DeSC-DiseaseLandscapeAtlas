# Energy Landscape Analysis of 51 Disease Families Reveals Comorbidity Archetypes and Transition Dynamics from 11 Million Patient Claims

Koki Tsuyuzaki^1,2,\*, Katsuyuki Chida^3,4, Tetsuo Ishikawa^3,5,6,7, and Eiryo Kawakami^3,4,5,8,\*

1: Data Science Core, Chiba University, 1-8-1 Inohana, Chuo, Chiba, Chiba 260-8670, JAPAN.
2: Omics AI Research Team, Advanced General Intelligence in Science Program (AGIS), TRIP Headquarters, RIKEN, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan.
3: Department of Artificial Intelligence Medicine, Graduate School of Medicine, Chiba University, 1-8-1 Inohana, Chuo, Chiba, Chiba 260-8670, JAPAN.
4: Predictive Medicine Special Project (PMSP), RIKEN Center for Integrative Medical Sciences (IMS), RIKEN, Yokohama, Kanagawa, JAPAN.
5: Division of Applied Mathematical Science, RIKEN Center for Interdisciplinary Theoretical and Mathematical Sciences (iTHEMS), RIKEN, Wako, Saitama, JAPAN.
6: Department of Extended Intelligence for Medicine, The Ishii-Ishibashi Laboratory, Keio University School of Medicine, 35 Shinanomachi, Shinjuku-ku, Tokyo, 160-8582, Japan.
7: Collective Intelligence Research Laboratory, Graduate School of Arts and Sciences, The University of Tokyo, 3-8-1 Komaba, Meguro-ku, Tokyo, 153-8902, Japan.
8: Department of Advanced Artificial Intelligence Medicine, Graduate School of Medicine, The University of Osaka, 2-2 Yamada-oka, Suita, Osaka 565-0871, Japan.

\*Corresponding authors: Koki Tsuyuzaki (koki.tsuyuzaki@gmail.com) and Eiryo Kawakami (eiryo.kawakami@chiba-u.jp)

---

## Summary

Comorbidity patterns shape clinical trajectories, yet systematic cross-disease comparisons remain rare. We apply energy landscape analysis — an Ising model framework from statistical mechanics — to administrative claims from approximately 11 million individuals in the Japanese DeSC database, constructing comorbidity landscapes for 51 disease families. For each family, a binary comorbidity matrix is reduced to seven dimensions by principal component analysis, discretized into 128 states, and analyzed to identify 1–4 stable basins per disease. In 27 of 51 families, lower energy is associated with higher mortality risk in at least half of basins. Among deceased patients in single-basin diseases, first and last basin assignments matched in a median of 98% of cases, indicating that basins function as persistent attractors. Cross-disease comparison of basin comorbidity codes recovers clinically coherent disease clusters without supervision. We classify the 51 families into six behavioral patterns (C1–C6), providing a descriptive atlas of comorbidity organization across diseases.

**Keywords:** energy landscape analysis, comorbidity, claims data, disease phenotyping, unsupervised learning, Ising model

---

## 1 Introduction

Administrative claims databases capture diagnoses, procedures, and prescriptions for millions of individuals over extended periods, providing an increasingly important resource for studying disease at the population level (Schneeweiss and Avorn, 2005). While individual comorbidities and their pairwise associations have been extensively catalogued (Barnett et al., 2012), the higher-order structure of comorbidity — how multiple conditions co-occur, form stable configurations, and evolve over time — has received less attention. Network approaches have begun to map disease–disease relationships at population scale (Hidalgo et al., 2009; Jensen et al., 2014), but methods that identify discrete comorbidity states and their stability remain underexplored. Understanding this structure matters because comorbidity burden is a major driver of healthcare utilization, mortality, and treatment complexity (Salisbury et al., 2011).

Energy landscape analysis (ELA) offers a framework for discovering stable states in high-dimensional binary data. Originally developed in statistical mechanics to describe spin systems and applied to protein folding (Bryngelson et al., 1995), the approach constructs an energy function over discrete states such that frequently occupied configurations receive low energy. Local minima of this energy surface define "basins of attraction" — groups of states separated by energy barriers. The maximum-entropy (Ising model) formulation, which infers pairwise couplings from observed state frequencies (Schneidman et al., 2006), has been applied to neural population activity (Watanabe et al., 2014; Ezaki et al., 2017) and gene expression states (Sasai and Wolynes, 2003), among other domains (for a tutorial, see Masuda et al., 2025).

Here, we present the Disease Landscape Atlas, a systematic application of ELA to 51 disease families defined by ICD-10 code prefixes in the DeSC Japanese claims database. By constructing separate landscapes for each disease family — excluding the index disease codes from the feature set — the resulting basins reflect disease-specific comorbidity patterns rather than population-level prevalence. The cohorts collectively span approximately 11 million individuals observed between April 2014 and February 2022, with patients potentially contributing to multiple disease families. For each family, we construct a binary comorbidity matrix from 6-month observation windows, apply PCA to reduce dimensionality to 7 components, binarize the scores, and estimate the energy landscape using the Ising model.

Our analysis addresses three questions. First, does the direction of the energy–outcome association hold broadly across diseases, and for which disease families does it fail? Second, do the comorbidity codes defining basins reflect known clinical relationships, and does cross-disease comparison of basin content recover established disease groupings? Third, how stable are basin assignments over time, and do different diseases exhibit distinct patterns of inter-basin transitions?

---

## 2 Results

### 2.1 Analysis design and disease-family-specific landscape construction

We defined 51 disease families by ICD-10 code prefixes (Figure 1; Methods), spanning cardiovascular, neoplastic, neuropsychiatric, respiratory, metabolic, musculoskeletal, and infectious categories. For each family, an index cohort was constructed by identifying all individuals with at least one claim matching the family's defining ICD-10 codes in the DeSC database. Disease family cohorts ranged from 47,332 individuals (inflammatory bowel disease) to 4,514,492 individuals (diabetes mellitus) (Table S1). Because disease families are defined by non-exclusive ICD-10 prefixes and individuals can carry multiple diagnoses, patients may contribute to more than one family; the cohort sizes therefore cannot be summed to obtain the total number of unique individuals.

For each cohort, we constructed a binary comorbidity matrix in which rows represent patient–observation window pairs (6-month windows advanced in 1-month steps, yielding up to 83 windows per patient) and columns represent ICD-10 codes, excluding the index disease codes themselves. PCA reduced this matrix to 7 principal components per disease, and scores were binarized (positive → +1, non-positive → −1), yielding 2^7 = 128 discrete spin states per disease family. The Ising model energy function was estimated from the observed state frequencies using the Landscaper software, and local energy minima were identified as basins of attraction. Across the 51 families, the number of basins ranged from 1 to 4 (17 families with 1 basin, 26 with 2, 6 with 3, and 2 with 4).

### 2.2 PCA captures age and sex gradients across disease families

Seven principal components captured a median of 12.5% of total variance across the 51 families (Figure 2c), yet they consistently captured clinically interpretable axes. PC1 correlated with patient age across nearly all disease families (median |r| = 0.41; Figure 2d), as illustrated by the diabetes cohort where older patients cluster at one pole of PC1 (Figure 2a). Sex separation was distributed across different PCs depending on the disease; in the urinary tract infection cohort, for example, PC7 separated female and male patients (|r| = 0.45; Figure 2b). The choice of 7 dimensions is a practical constraint: ELA requires enumerating all 2^dim states, making 7 the maximum (128 states) for the Landscaper software. The low cumulative variance is expected given the high dimensionality of ICD-10 code space; full covariance analyses of similar claims matrices suggest that hundreds of components are needed to account for the majority of variance. The 7 components nonetheless preserve the major clinical gradients used for landscape construction. The remaining components (typically PC2–PC6) captured disease-specific comorbidity patterns whose loadings varied across families; systematic interpretation of these axes is beyond the scope of this study.

### 2.3 Energy is associated with mortality and hospitalization across disease families

To assess whether the energy landscape captures clinically relevant structure, we tested the association between within-basin standardized energy and all-cause mortality using Cox proportional hazards models (adjusted for age and sex) within each basin separately (Figure 3; Methods §4.5). Mortality analyses were restricted to patients enrolled in National Health Insurance (NHI) or Late-Stage Elderly insurance, which record death events (approximately 40–50% of each cohort; see Methods). A disease family was classified as satisfying the per-basin majority criterion if HR < 1 (lower energy associated with higher hazard) held in at least half of its basins. In 27 of 51 disease families, this criterion was satisfied for mortality. An analogous analysis using disease-specific hospitalization as the endpoint, which does not require death flag availability and thus used the full population, yielded a consistent direction in 35 of 51 families (Figure S1).

The twenty-four disease families that did not satisfy the mortality criterion were acute renal failure (2 basins, 0/2 HR < 1), anxiety (2 basins, 0/2), autoimmune/connective tissue disease (1 basin, 0/1), benign prostatic hyperplasia (3 basins, 1/3), cirrhosis/fibrosis (1 basin, 0/1), COPD (1 basin, 0/1), delirium (2 basins, 0/2), depression (2 basins, 0/2), dyslipidemia (2 basins, 0/2), epilepsy (2 basins, 0/2), gout (1 basin, 0/1), hip fracture (1 basin, 0/1), liver cancer (3 basins, 0/3), pancreatic cancer (1 basin, 0/1), Parkinson disease (3 basins, 1/3), pneumonia (2 basins, 0/2), prostate cancer (2 basins, 0/2), rheumatoid arthritis (2 basins, 0/2), schizophrenia (2 basins, 0/2), sepsis (3 basins, 0/3), thyroid disorder (4 basins, 0/4), tuberculosis (1 basin, 0/1), urinary tract infection (3 basins, 1/3), and viral hepatitis (2 basins, 0/2). These diseases include non-terminal chronic conditions with low mortality rates (anxiety, dyslipidemia, gout, benign prostatic hyperplasia, prostate cancer), neuropsychiatric conditions (depression, schizophrenia, epilepsy, delirium), diseases with acute presentation or high lethality (sepsis, acute renal failure, urinary tract infection), and conditions with progressive or organ-specific courses (cirrhosis/fibrosis, hip fracture, pancreatic cancer, Parkinson disease, pneumonia, viral hepatitis). The non-terminal diseases may produce comorbidity landscapes where energy reflects healthcare utilization or diagnostic complexity rather than mortality risk; the acute diseases may lack sufficient longitudinal comorbidity information for severity-aligned landscapes (Figure 3).

We emphasize that this criterion is a directional consistency check, not a predictive performance evaluation. The large sample sizes (up to 4.5 million patients per family) mean that even small and clinically negligible effect sizes can reach statistical significance; the hazard ratios should be interpreted as indicating the direction and approximate magnitude of the energy–outcome association, not as evidence of predictive utility.

### 2.4 Basin comorbidity codes recover clinically coherent disease clusters

To assess whether ELA basins capture clinically meaningful comorbidity structure, we compared the ICD-10 code compositions of basins across disease families. For each basin, we selected the 30 most frequent ICD-10 codes by within-basin observation-window frequency. For each pair of diseases, we computed the maximum Jaccard similarity between any two of their respective basins' code sets. Schizophrenia was excluded from this analysis because its antipsychotic-dominated code profile was an outlier in the similarity space, yielding 50 disease families (Figure 4; Methods).

Hierarchical clustering of the resulting distance matrix (1 − maximum Jaccard; average linkage with optimal leaf ordering) revealed two major groups. An Acute/Frailty cluster contained 3 diseases: acute renal failure, aspiration pneumonia, and sepsis. The remaining 47 diseases formed a broad Chronic/Lifestyle cluster. Within the Chronic/Lifestyle cluster, clinically coherent sub-clusters emerged without supervision: neurodegenerative diseases (dementia, Parkinson disease, epilepsy), a frailty-adjacent group (delirium, pneumonia, malnutrition, hip fracture) positioned nearest the Acute/Frailty cluster, prostate diseases (benign prostatic hyperplasia, prostate cancer), pulmonary fibrosis diseases (interstitial lung disease, lung cancer), musculoskeletal diseases (osteoarthritis, osteoporosis), and a hepatic group (cirrhosis/fibrosis, liver cancer), among others.

These groupings, derived solely from the frequency of comorbidity codes within ELA basins, recapitulate known disease relationships without supervision. Certain comorbidity codes — such as chronic periodontitis and hyperlipidemia — discriminated basins across many diseases, suggesting shared comorbidity axes (Figure S4). Within-disease basin similarity varied widely across diseases (Figure S3). For example, neurodegenerative diseases showed high within-disease Jaccard values (basins sharing similar comorbidity profiles), whereas diseases like diabetes showed low values, indicating that its basins harbor distinct comorbidity profiles — in some cases more similar to basins of other diseases than to each other.

### 2.5 Longitudinal basin transitions and six behavioral patterns (C1–C6)

To characterize the temporal stability of basin assignments, we tracked the basin membership of deceased patients across observation windows preceding death (Figure 5). Three observations emerged.

First, basins function as persistent attractors. Among deceased patients whose first observed state fell within a basin, a median of 98% had the same basin assignment at their first and last observation in single-basin diseases (17 families; Figure 5b). This first-to-last match rate permits intermediate departures and returns — it does not require continuous residence in the same basin. In diseases with 3–4 basins, the first-to-last match rate was approximately 69%, with approximately 30% moving to a different basin by the final observation. Exit from basins to non-basin states was rare (median below 3% across basin-count groups). Second, non-basin states are transient: among deceased patients starting outside a basin, a median of 74–90% ended in a basin (Figure 5c). Third, the fraction of patients remaining entirely outside basins throughout follow-up (Stable Non-Basin in Figure 5a) was generally below 3% (median 0.4%, maximum 4.4% in obesity). These statistics are medians across diseases within each basin-count group, computed from deceased patients only.

To organize the diversity of landscape behaviors across the 51 families, we defined a descriptive classification based on three axes: number of basins, direction of the energy–mortality association (per-basin majority criterion), and transition asymmetry — the degree to which inter-basin transition rates differ between directions (e.g., Basin 1→2 vs. 2→1; see Methods) — computed from deceased patients (Table 1; Figure S2). We label the six resulting categories C1–C6 for conciseness; the numbering is arbitrary and does not imply a severity or progression ordering. This classification is descriptive rather than algorithmic: it was defined post hoc to summarize observed patterns, not discovered automatically or validated on independent data.

**C1** (10 families): Single basin, majority criterion satisfied. These diseases — including stroke, ischemic heart disease, colorectal cancer, lung cancer, and metastatic cancer — have a single stable comorbidity configuration with a clear severity gradient within it.

**C2** (7 families): Single basin, majority criterion not satisfied. COPD, autoimmune/connective tissue disease, tuberculosis, cirrhosis/fibrosis, gout, hip fracture, and pancreatic cancer have a single basin in which lower energy is not associated with worse outcomes.

**C3** (2 families): Multiple basins, majority criterion satisfied, low transition asymmetry (<0.1). CKD and dementia: basins exchange patients at roughly equal rates. In dementia, for instance, the full-period per-step transition rates are Basin 1→2 = 5.7% and Basin 2→1 = 5.9% (asymmetry = 0.02; 241,849 and 292,895 transitions respectively).

**C4** (6 families): Multiple basins, majority criterion satisfied, moderate transition asymmetry (0.1–0.45). Asthma, atrial fibrillation, heart failure, hypertension, malnutrition, and sleep disorders exhibit asymmetric conditional transition probabilities between basins.

**C5** (9 families): Multiple basins, majority criterion satisfied, high asymmetry (≥0.45). Anemia, malignancy, diabetes, and osteoarthritis typify this pattern: inter-basin transitions are strongly directional when they occur.

**C6** (17 families): Multiple basins, majority criterion not satisfied. Anxiety, benign prostatic hyperplasia, delirium, depression, dyslipidemia, epilepsy, liver cancer, Parkinson disease, pneumonia, prostate cancer, rheumatoid arthritis, schizophrenia, sepsis, thyroid disorder, acute renal failure, urinary tract infection, and viral hepatitis exhibit structured basin organization but the expected energy–outcome direction does not hold in the majority of basins. As a representative C6 example, sepsis showed three basins with divergent survival curves for both all-cause mortality (Figure S5) and disease-specific hospitalization (Figure S6), and a structured state occupancy pattern before death (Figure S7).

The asymmetry boundaries (<0.1, 0.1–0.45, ≥0.45) were determined from a gap analysis of the corrected asymmetry distribution (Figure S2). Non-basin states were excluded from asymmetry computation to prevent inflation of basin counts and transition rates.

### 2.6 Diabetes: C5 pattern with cardiovascular and dental basins

The diabetes cohort (n = 4,514,492; ICD-10 E10–E14) yielded four basins with a C5 pattern (Figure 6). The disconnectivity graph — which shows the energy barriers separating basins, with deeper merges indicating more closely related basins — reveals two pairs: Basin 1 (E = −1.01) and Basin 3 (E = −0.95) merge at the lowest connection energy, forming the most stable pair, while Basin 2 (E = −0.21) and Basin 4 (E = −0.35) are shallower and merge at higher barriers (Figure 6c). This hierarchical structure reflects the landscape topology but does not itself indicate clinical severity.

Basin 1 (n = 1,926,087) was characterized by cardiovascular comorbidity codes ranked by log-odds difference score (Methods): pure hypercholesterolemia (E780; 26.6% of observation windows), essential hypertension (I10; 37.1%), old myocardial infarction (I252; 1.4%), unstable angina (I200; 1.0%), and angina pectoris (I209; 10.2%) (Figure 6a). Basin 2 (n = 186,913) was characterized by dental codes: dental caries (K029; 17.9%), surgical complications (T888; 21.5%), loss of teeth (K081; 12.2%), and periodontitis-related codes. Among the mortality-analysis population (NHI and Late-Stage Elderly enrollees; n = 1,839,848), Basin 4 had the lowest death rate (6.1%; 7,069/115,104), while Basin 2 had the highest (12.9%; 11,323/88,051); Basin 1 had 9.1% (59,983/658,494) and Basin 3 had 10.4% (101,795/978,199) (Figure 6b). The per-basin majority criterion was satisfied (2 of 4 basins with HR < 1): Basin 1 (HR = 0.844) and Basin 4 (HR = 0.967) showed the expected direction, while Basin 2 (HR = 1.193) and Basin 3 (HR = 1.045) were reversed. Temporal basin composition among deceased patients was largely stable across observation windows preceding death (Figure 6d).

The transition asymmetry was 0.56, placing diabetes among the C5 diseases. In the transition matrix among deceased patients (last 20 observation positions), Conditional transition probabilities were strongly asymmetric: Basin 2→Basin 1 (11.4%) vs. Basin 1→Basin 2 (1.9%), and Basin 4→Basin 3 (13.7%) vs. Basin 3→Basin 4 (0.7%). However, because Basin 1 and Basin 3 are much larger, the absolute transition counts were roughly comparable (e.g., Basin 2→1: 22,923 vs. Basin 1→2: 21,733; Basin 4→3: 14,544 vs. Basin 3→4: 12,759). The asymmetry in conditional probabilities therefore reflects basin size differences rather than a large net patient flow. The co-occurrence patterns describe associations, not causal relationships.

### 2.7 Dementia: C3 pattern with divergent mortality between basins

The dementia cohort (n = 933,285; ICD-10 F00–F03, G30) yielded two basins with a C3 pattern. Over the full observation period, Basin 1→2 transitions occurred at 5.7% per step and Basin 2→1 at 5.9% per step (asymmetry = 0.02; 241,849 and 292,895 transitions respectively), indicating balanced exchange between basins despite their markedly different mortality profiles.

Among the mortality-analysis population (NHI and Late-Stage Elderly enrollees; n = 449,369), the two basins showed sharply different mortality: Basin 1 (n = 183,908) had a death rate of 33.6%, compared with 17.2% in Basin 2 (n = 265,461). The energy–mortality association held in Basin 1 (HR = 0.855, p < 0.001) but was reversed in Basin 2 (HR = 1.132, p < 0.001). Dementia satisfies the per-basin majority criterion (1 of 2 basins HR < 1) but has a split energy–outcome pattern.

The two basins showed different comorbidity profiles. Basin 1 (n = 370,219 in the Cox-eligible cohort; mean age 83.0) was characterized by constipation (79%), xerosis cutis (39%), dermatitis (38%), and heart failure (39%). Basin 2 (n = 511,818; mean age 80.8) was characterized by hypertension (78%), diabetes (48%), and hyperlipidemia (48%). Basin 1's frailty-associated profile aligns with its higher mortality, while Basin 2's metabolic profile (hypertension, diabetes, hyperlipidemia) represents a less acutely life-threatening comorbidity pattern. Despite this mortality gap, inter-basin exchange over the full period is nearly balanced (5.7% vs. 5.9% per step; asymmetry = 0.02), suggesting coexisting subtypes rather than a progressive trajectory.

The step-by-step transition matrices (Figure 7) illustrate how transition dynamics vary across the C3–C5 categories. In dementia (C3; Figure 7a), deceased-patient transition rates computed from the last 20 observation positions before death are Basin 1→2 = 3.9% and Basin 2→1 = 10.3% — more asymmetric than the full-period values (5.7% and 5.9%) because the time window is restricted to the final observation positions. In heart failure (C4; Figure 7b), Basin 3→1 = 11.1% vs. Basin 1→3 = 2.9%, reflecting moderate asymmetry. In diabetes (C5; Figure 7c), transition asymmetry is pronounced across all four basins (e.g., Basin 2→1 = 11.4% vs. Basin 1→2 = 1.9%). The C3–C5 classification is based on full-period asymmetry values; Figure 7 shows the corresponding patterns in the pre-death window.

### 2.8 Tuberculosis: a structurally simple C2 landscape

The TB cohort yielded a single basin, placing it in the C2 category. Among TB patients in the mortality-analysis population (n = 148,146), the energy–mortality association was reversed (HR = 1.143; 95% CI 1.130–1.157), though the mean energy difference between deceased and surviving patients was small (−0.080 vs. −0.091) (Figure 8). TB patients had among the fewest observation windows per patient in the atlas (9.8, compared with a median of ~20 across all 51 families), reflecting short insurance observation periods.

Among deceased TB patients, 82.6% remained in Basin 1 across all observed windows, while 1.5% never entered Basin 1. Figure 8a shows the fraction of patients in each state category at each time point before the last retained observation.

The landscape was partitioned by a single PC axis: PC6 had value −1 in the basin state and +1 in the top-10 most-frequent non-basin states, while PCs 1–5 and 7 varied freely (Figure 8b). The ICD-10 codes contributing most to this discriminating axis are shown in Figure 8c. In the Jaccard clustering (Figure 4), TB grouped with respiratory diseases (interstitial lung disease, lung cancer, COPD) rather than with the Acute/Frailty cluster, consistent with its comorbidity profile reflecting chronic pulmonary conditions. The short observation periods likely limit the temporal comorbidity variation available for landscape construction, though other diseases with similarly few observation windows (e.g., aspiration pneumonia, 9.5 windows/patient) do satisfy the majority criterion.

---

## 3 Discussion

We have presented the Disease Landscape Atlas, a systematic application of energy landscape analysis to 51 disease families in a Japanese claims database of approximately 11 million individuals. The atlas reveals that ELA basins capture clinically interpretable comorbidity structure in the majority of disease families and that this structure can be organized into six behavioral patterns.

**Basins as persistent comorbidity configurations.** Basin assignments show high temporal persistence among deceased patients: a median of 98% first-to-last concordance in single-basin diseases and 69% in diseases with 3–4 basins (Figure 5b), while non-basin states are transient (74–90% of patients starting outside basins ending in one). The degree of inter-basin exchange ranges from the near-balanced exchange in dementia (C3; asymmetry = 0.02) to the highly asymmetric transitions of anemia (C5; ~110:1 ratio), consistent with clinical expectations — in anemia, for example, the asymmetry aligns with progressive renal decline driving chronic anemia through erythropoietin deficiency (Babitt and Lin, 2012).

**Cross-disease validity of basin structure.** The Jaccard clustering (Figure 4) recovered clinically coherent disease groups — neurodegenerative diseases, prostate diseases, musculoskeletal diseases, and a tight acute/frailty core — without supervision, echoing groupings identified by network approaches using different data and methods (Hidalgo et al., 2009; Jensen et al., 2014). Frailty-associated conditions (delirium, pneumonia, malnutrition, hip fracture) fell at the boundary between the Chronic/Lifestyle and Acute/Frailty clusters, consistent with their dual character.

**Energy–outcome association.** The per-basin majority criterion was satisfied in 27 of 51 families for mortality and 35 of 51 for hospitalization. The 24 families not satisfying the mortality criterion span non-terminal chronic conditions (anxiety, gout, benign prostatic hyperplasia, prostate cancer), neuropsychiatric conditions (depression, schizophrenia, epilepsy, delirium), and diseases with short observation periods or high lethality (sepsis, acute renal failure, tuberculosis). The dementia case study illustrates that basins with markedly different mortality (33.6% vs. 17.2%) can exchange patients at nearly equal rates (asymmetry = 0.02), suggesting coexisting subtypes rather than progressive stages.

**The six-pattern classification.** The C1–C6 categories provide a descriptive vocabulary for the behavioral diversity of comorbidity landscapes. The classification was defined post hoc to summarize observed variation; the numbering is arbitrary and does not imply severity ordering, and the asymmetry boundaries were chosen from distributional gaps rather than external criteria.

**Limitations.** ICD-10 codes from administrative claims are recorded for billing purposes and may not reflect clinical reality (O'Malley et al., 2005). The overlapping observation windows (5 of 6 months shared) may inflate apparent temporal stability. The 7-dimensional PCA reduction captures only 12.5% of total variance, and the 128 discrete states are a coarse approximation of the comorbidity space. The Cox analysis uses each patient's first eligible observation window as the representative record (see Methods), which may underestimate the energy–mortality association for diseases where comorbidity burden accumulates slowly. Proportional hazards assumption testing, competing risk analysis, and multiple testing correction were not performed; we focus on effect direction rather than statistical significance. Trajectory and transition analyses are restricted to deceased patients. Mortality analyses are restricted to NHI and Late-Stage Elderly enrollees (approximately 40–50% of each cohort), and results may not generalize to employer-based insurance enrollees. The mortality and hospitalization analyses use different population subsets (e.g., dementia: n = 449,369 vs. n = 858,097), so absolute rates between endpoints should be compared with caution.

**Outlook.** The Disease Landscape Atlas demonstrates that ELA can extract stable, interpretable, and cross-disease-comparable comorbidity structures from claims data at scale. Potential applications include identifying patient subgroups for stratified analysis, generating hypotheses about comorbidity-driven disease trajectories, and comparing comorbidity architecture across populations or healthcare systems. Because a patient's current comorbidity codes can be projected onto the stored eigenvectors to obtain a basin assignment and energy value, and the inter-basin transition matrices can serve as a Markov model, the atlas could in principle be extended to estimate individual-level trajectory and outcome probabilities. However, the transition matrices in this study were estimated from deceased patients only; whether they generalize to the full population remains to be tested, and predictive accuracy would require formal evaluation. Validation in independent cohorts and formal assessment of the relationship between energy barriers and transition dynamics represent important next steps.

---

## 4 Methods

### 4.1 Data source and construction of disease-specific matrices

The analysis used administrative claims data from the DeSC database, a Japanese claims and health checkup database, covering the period from April 2014 to February 2022. The database contains diagnosis codes (ICD-10, 2013 Japanese modification), demographic information (birth year-month, sex, insurer type), observation period (start and end of coverage), death flag (available for National Health Insurance and Late-Stage Elderly enrollees; see Methods §4.5), and hospitalization records.

We defined 51 disease families based on ICD-10 code prefixes (Supplementary Table S1). Each family was specified by one or more ICD-10 prefixes; for example, the dementia family included F00–F03 and G30, while the diabetes family included E10–E14. For each disease family, an index cohort was defined as all individuals with at least one claims record matching the family's ICD-10 prefix codes during the observation period. Observation windows of 6 months' duration were constructed with a sliding step of 1 month, yielding up to 83 windows per patient across the full observation period.

A binary comorbidity matrix was constructed for each cohort, with rows representing patient–observation window pairs and columns representing ICD-10 diagnosis codes. The matrix entry was 1 if the code appeared in at least one claim within that observation window and 0 otherwise. The ICD-10 codes defining the index disease were excluded from the column set for each family, so that the resulting landscape reflects comorbidity patterns rather than the index disease itself. After construction, the matrix was stored in coordinate (COO) format with zstd compression.

### 4.2 Principal component analysis

For each disease family, the comorbidity matrix was subjected to exact, out-of-core PCA using OnlinePCA.jl (Tsuyuzaki, 2026; Tsuyuzaki et al., 2020) with a chunk size of 300,000 rows. The first 7 principal components were retained for all families. The choice of 7 dimensions was driven by the Landscaper software's requirement to enumerate all 2^dim states; 7 dimensions yield 128 states. The median cumulative variance explained by 7 PCs was 12.5% (Figure 2c; computed from patient-representative scores; eigenvalue-based calculation yields 13.3%).

### 4.3 Binarization of principal component scores

Each PCA score was binarized using a zero threshold: values greater than 0 were assigned B = 1, and values less than or equal to 0 were assigned B = 0. The corresponding Ising spin variable was defined as S = 2B − 1, yielding S ∈ {−1, +1} for each of the 7 components. Each patient–observation window was thereby assigned one of 2^7 = 128 discrete states. The Landscaper software expects a sparse binary matrix where non-zero entries indicate S = +1 and absent entries indicate S = −1.

### 4.4 Energy landscape analysis using Landscaper

Energy landscape estimation was performed using Landscaper (https://github.com/chiba-ai-med/Landscaper), which fits an Ising model to the observed state frequencies and computes the energy E(σ) for each of the 128 states σ. States with lower energy correspond to higher observation probability. Local energy minima were identified, and each state was assigned to the basin of the nearest local minimum under steepest-descent dynamics. The SubGraph.tsv output records basin labels for states that belong to a non-trivial basin. In 47 of 51 diseases, SubGraph.tsv excludes isolated local minima (states with no lower-energy neighbor in the StatusNetwork), yielding fewer than 128 rows; in the remaining 4 diseases (dementia, delirium, diabetes, sepsis), all 128 states are included. Basin.tsv records the state index of each local minimum. The disconnectivity graph shows the connection energies at which basin pairs merge. Landscaper was run with a single thread per disease.

### 4.5 Basin assignment and clinical outcome analysis

**Representative record selection.** Because each patient contributed multiple observation windows, a single representative record per patient was required for Cox analysis. For each patient, the first observation window satisfying all of the following criteria was selected: (1) the full 6-month window falls within the patient's insurance observation period, (2) the window falls within the data coverage period (April 2014 to February 2022), (3) the window contains at least one claim matching the index disease, and (4) the patient's state has a basin assignment. This record's state, basin assignment, and energy value served as the patient-level summary. Follow-up time begins at the end of this representative window, ensuring that no survival time accrues before the comorbidity state is measured.

**Cox proportional hazards models.** For each disease family and each basin, a Cox proportional hazards model was fitted:

> h(t) = h₀(t) × exp(β₁ E_z + β₂ age + β₃ sex)

where E_z is the energy value standardized (Z-score) within each basin, age is calculated from birth year-month to the representative window start, and sex is coded as a factor.

*Mortality.* The time variable is the interval from the end of the representative window to the patient's insurance coverage end (converted to years). The DeSC database records a death flag for patients enrolled in National Health Insurance (NHI) and Late-Stage Elderly insurance; it is not recorded for employer-based insurance enrollees. All mortality analyses were restricted to patients for whom the death flag was available. Patients who died before or during the representative window, or whose coverage ended before it, were excluded. Under the Japanese NHI and Late-Stage Elderly systems, insurance coverage terminates upon death, so the coverage end month approximates the month of death for these enrollees.

*Disease-specific hospitalization.* The time origin is likewise the end of the representative window. Patients whose first qualifying admission occurred at or before this time were excluded. The event is first admission where the DPC primary diagnosis matches the index disease family's ICD-10 codes, occurring after the representative window end and within the data period (April 2014 to February 2022). Death is not explicitly modeled as a competing risk. The hospitalization analysis uses all eligible patients (not restricted to death-flag-available enrollees) but excludes those with a qualifying admission before the representative window end, yielding different sample sizes than the mortality analysis (e.g., dementia: n = 449,369 for mortality vs. n = 858,097 for hospitalization).

**Per-basin majority criterion.** A disease family was classified as satisfying the criterion if HR < 1 held in at least half of its basins. This includes ties: for a disease with 2 basins, 1 of 2 suffices; for 4 basins, 2 of 4 suffices.

**Kaplan–Meier curves.** Kaplan–Meier survival curves were constructed for inter-basin comparison, stratified by basin assignment from the representative record. For visualization, energy tertile groups (low, medium, high) within each basin were used to illustrate intra-basin variation.

### 4.6 Comorbidity profiles and cross-disease similarity

**Log-odds difference scoring.** To characterize each basin's comorbidity profile, we computed a log-odds difference score for each ICD-10 code j in basin b:

> score(b, j) = logit( (c_bj + α) / (M_b + 2α) ) − logit( (c_all_j + α) / (M_all + 2α) )

where c_bj is the number of observation windows in basin b containing code j, M_b is the total count of all code occurrences across observation windows in basin b, c_all_j and M_all are the corresponding totals across all basins (including basin b itself), α = 1 is a smoothing parameter, and logit(p) = ln(p / (1 − p)). Codes were filtered to those with within-basin observation-window frequency > 1% and at least 100 occurrences.

**Jaccard similarity and cross-disease clustering.** For each basin, the top 30 ICD-10 codes were selected by within-basin observation-window frequency (not by log-odds difference score). For each pair of disease families, the maximum Jaccard index between any two of their basins' top-30 code sets was computed. The resulting similarity matrix was converted to a distance matrix (1 − Jaccard) and subjected to average-linkage hierarchical clustering with optimal leaf ordering. Schizophrenia was excluded, leaving 50 disease families.

### 4.7 Longitudinal trajectory and transition analysis

**Trajectory-to-death analysis.** For each deceased patient in each disease family, we tracked the patient's state, energy, and basin assignment across all retained observation windows, indexed by ordinal position counting backward from the last retained observation (position 0 at the final retained window). The ordinal positions reflect the sequence of available records, not calendar-time intervals. The full trajectory (up to 83 positions per patient) was stored and used for the full-period transition rate computation underlying the C3–C5 classification (Table 1, Figure S2). For visualization (trajectory heatmaps, basin composition plots) and the transition comparison in Figure 7, only the last 20 observation positions were used.

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

6. Jensen AB, Moseley PL, Oprea TI, et al. Temporal disease trajectories condensed from population-wide registry data covering 6.2 million patients. *Nature Communications*. 2014;5:4022. doi:10.1038/ncomms5022. PMID: 24959948.

7. Masuda N, Islam S, Aung ST, Watanabe T. Energy landscape analysis based on the Ising model: tutorial review. *PLOS Complex Systems*. 2025;2(5):e0000039. doi:10.1371/journal.pcsy.0000039.

8. Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33. doi:10.12688/f1000research.29032.2.

9. O'Malley KJ, Cook KF, Price MD, Wildes KR, Hurdle JF, Ashton CM. Measuring diagnoses: ICD code accuracy. *Health Services Research*. 2005;40(5 Pt 2):1620–1639. doi:10.1111/j.1475-6773.2005.00444.x. PMID: 16178999.

10. Salisbury C, Johnson L, Purdy S, Valderas JM, Montgomery AA. Epidemiology and impact of multimorbidity in primary care: a retrospective cohort study. *British Journal of General Practice*. 2011;61(582):e12–e21. doi:10.3399/bjgp11X548929.

11. Sasai M, Wolynes PG. Stochastic gene expression as a many-body problem. *Proceedings of the National Academy of Sciences*. 2003;100(5):2374–2379. doi:10.1073/pnas.2627987100.

12. Schneeweiss S, Avorn J. A review of uses of health care utilization databases for epidemiologic research on therapeutics. *Journal of Clinical Epidemiology*. 2005;58(4):323–337. doi:10.1016/j.jclinepi.2004.10.012. PMID: 15862718.

13. Schneidman E, Berry MJ II, Segev R, Bialek W. Weak pairwise correlations imply strongly correlated network states in a neural population. *Nature*. 2006;440(7087):1007–1012. doi:10.1038/nature04701. PMID: 16625187.

14. Tsuyuzaki K. OnlinePCA.jl: A Julia Package for Out-of-core and Sparse Principal Component Analysis. *Journal of Open Source Software*. 2026;11(117):9343. doi:10.21105/joss.09343.

15. Tsuyuzaki K, Sato H, Sato K, Nikaido I. Benchmarking principal component analysis for large-scale single-cell RNA-sequencing. *Genome Biology*. 2020;21:9. doi:10.1186/s13059-019-1900-3. PMID: 31955711.

16. Watanabe T, Masuda N, Megumi F, Kanai R, Rees G. Energy landscape and dynamics of brain activity during human bistable perception. *Nature Communications*. 2014;5:4765. doi:10.1038/ncomms5765. PMID: 25163855.

---

## Figure legends

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 1. Overview of the Disease Landscape Atlas analysis workflow.

(a) From the DeSC claims database, 51 disease-family cohorts are defined by ICD-10 index codes. For each cohort, a binary comorbidity matrix is constructed in which rows represent patient–observation window pairs and columns represent ICD-10 codes excluding the index disease. PCA reduces the matrix to 7 dimensions, scores are binarized to ±1, yielding 2^7 = 128 discrete states. Energy landscape analysis assigns each state an energy value and groups states into basins of attraction (1–4 per disease). (b) Within each basin, Cox proportional hazards models test whether lower energy is associated with higher mortality or hospitalization risk (HR < 1 per 1 SD increase in energy, adjusted for age and sex). The per-basin majority criterion (HR < 1 in at least half of basins) is satisfied in 27 of 51 families for mortality and 35 of 51 for hospitalization.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 2. PCA captures age and sex axes across disease families.

(a) Diabetes cohort projected onto PC1 vs. PC2, colored by age. PC1 correlates with age across most diseases (median |r| = 0.41). (b) Urinary tract infection cohort projected onto PC1 vs. PC7, colored by sex. PC7 separates female and male patients (|r| = 0.45). (c) Cumulative variance explained by the first 7 PCs for all 51 disease families (grey) and their median (red). Seven PCs capture a median of 12.5% of total variance. (d) Pearson correlation of each PC with age (left) and sex (right) across 51 disease families, sorted by PC1–age correlation. PC1 universally captures an age gradient; sex associations are distributed across different PCs per disease.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 3. Atlas of 51 disease families by observation density, death rate, and energy–mortality association.

Each point represents one disease family. The x-axis shows the mean number of patient–observation window rows per patient (observation density), and the y-axis shows the cohort death rate during follow-up. Point color indicates whether energy is associated with mortality in the expected direction (HR < 1) in at least half of basins (blue) or not (red). The per-basin majority criterion is satisfied in 27 of 51 disease families.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 4. Disease clustering by basin code similarity recovers clinically valid disease groups.

Hierarchical clustering of 50 disease families (schizophrenia excluded) by basin code similarity. For each pair of diseases, the maximum Jaccard similarity between any two basins — using each basin's top-30 most frequent ICD-10 codes by within-basin observation-window frequency — is computed. Average-linkage clustering with optimal leaf ordering reveals two major groups: an Acute/Frailty cluster (3 diseases: acute renal failure, aspiration pneumonia, sepsis) and a broad Chronic/Lifestyle cluster (47 diseases). Within the Chronic/Lifestyle cluster, clinically coherent sub-clusters emerge without supervision, including a frailty-adjacent group at the boundary nearest the Acute/Frailty cluster.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 5. Basin assignments persist among deceased patients, with basin-count-dependent transition rates.

(a) Distribution of deceased patients' trajectory categories by number of basins. Each bar shows the median across diseases in that group. In 1-basin diseases (17 families), a median of 87% remain in the same basin throughout follow-up; in 2+ basin diseases, roughly half undergo state transitions. Stable Non-Basin is generally below 3% (median 0.4%). (b) Among deceased patients starting in a basin, a median of 98% had the same basin at first and last observation in 1-basin diseases; this drops to approximately 69% in 3–4 basin diseases, where approximately 30% move to a different basin. This is an endpoint comparison that permits intermediate departures and returns. (c) Among deceased patients starting outside a basin, a median of 74–90% ended in a basin. n = 24 for the 2-basin group, as 2 of 26 two-basin families had no non-basin-starting deceased patients.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 6. Diabetes energy landscape with cardiovascular and dental basins (C5).

(a) Top comorbidity codes in Basin 1 and Basin 2, ranked by log-odds difference score. Bar color indicates within-basin observation-window frequency. Patient counts are from the full cohort: Basin 1 (n = 1,926,087) is characterized by cardiovascular codes: pure hypercholesterolemia (E780; 26.6%), essential hypertension (I10; 37.1%), old myocardial infarction (I252; 1.4%), unstable angina (I200; 1.0%), and angina pectoris (I209; 10.2%). Basin 2 (n = 186,913) is characterized by dental codes: dental caries (K029; 17.9%), surgical complications (T888; 21.5%), and loss of teeth (K081; 12.2%). (b) Kaplan–Meier survival curves for all four basins, restricted to the mortality-analysis population (NHI and Late-Stage Elderly enrollees with death-flag availability; n = 1,839,848). Death rates (deaths/N per basin): Basin 4 = 6.1% (7,069/115,104), Basin 1 = 9.1% (59,983/658,494), Basin 3 = 10.4% (101,795/978,199), Basin 2 = 12.9% (11,323/88,051). The per-basin majority criterion is satisfied (2 of 4 basins with HR < 1): Basin 1 (HR = 0.844) and Basin 4 (HR = 0.967) show the expected direction, while Basin 2 (HR = 1.193) and Basin 3 (HR = 1.045) are reversed. (c) Disconnectivity graph. Basin 1 (E = −1.01) and Basin 3 (E = −0.95) merge at the lowest connection energy, forming the deepest pair; Basin 2 (E = −0.21) and Basin 4 (E = −0.35) join at higher barriers. (d) Temporal basin composition among deceased patients at successive observation positions before the last retained record.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 7. Basin transition probabilities before death in dementia, heart failure, and diabetes.

Step-by-step basin transition matrices computed from the last 20 observation positions of deceased patients for three representative diseases: (a) dementia (C3, 2 basins), (b) heart failure (C4, 2 basins), and (c) diabetes (C5, 4 basins). Rows represent the origin basin ("From") and columns represent the destination basin ("To"); each cell shows the conditional probability of transitioning from the row basin to the column basin at each observation step. Diagonal cells indicate retention within the same basin. Dementia (a): Basin 1→2 = 3.9%, Basin 2→1 = 10.3%. Heart failure (b): Basin 1→3 = 2.9%, Basin 3→1 = 11.1%. Diabetes (c): Basin 1→2 = 1.9%, Basin 2→1 = 11.4%. These transition rates differ from the full-period values used for the C3–C5 classification (Section 2.5, Table 1) because the time window differs: the figure uses the last 20 observation positions before death, whereas the full-period calculation uses all observation windows. Both analyses are restricted to deceased patients.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure 8. Tuberculosis — ELA limitations for low-density single-basin diseases (C2).

(a) Temporal basin composition among deceased TB patients at successive observation positions. (b) PC sign pattern comparison between the basin state and the top-10 most-frequent non-basin states. PC6 is the sole axis that consistently differs: PC6 = −1 in the basin state vs. PC6 = +1 in all 10 examined non-basin states. This comparison is limited to the selected states. (c) PC6 loading bar chart showing the ICD-10 codes contributing to the discriminating axis. Among the mortality-analysis population (n = 148,146), Cox analysis yields HR = 1.143 (95% CI: 1.130–1.157), indicating reversed energy–mortality direction.

---

```{=openxml}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="11906" w:h="16838" w:orient="portrait"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

## Table 1. Six-pattern classification of 51 disease families

| Category | Basins | Cox | Asymmetry | n | Diseases |
|----------|--------|-----|-----------|---|----------|
| C1 | 1 | Satisfied | — | 10 | Breast ca., Colorectal ca., ILD, IHD, Liver disease, Lung ca., Metastatic ca., Obesity, Stomach ca., Stroke |
| C2 | 1 | Not satisfied | — | 7 | Autoimmune conn. tissue, Cirrhosis, COPD, Gout, Hip fracture, Pancreatic ca., TB |
| C3 | Multi | Satisfied | <0.1 | 2 | CKD, Dementia |
| C4 | Multi | Satisfied | 0.1–0.45 | 6 | Asthma, AF, Heart failure, HTN, Malnutrition, Sleep disorders |
| C5 | Multi | Satisfied | ≥0.45 | 9 | Anemia, Asp. pneumonia, Diabetes, Hem. malignancy, IBD, Malignancy, OA, Osteoporosis, Peptic ulcer |
| C6 | Multi | Not satisfied | Any | 17 | AKI, Anxiety, BPH, Delirium, Depression, Dyslipidemia, Epilepsy, Liver ca., Parkinson, Pneumonia, Prostate ca., RA, Schizophrenia, Sepsis, Thyroid, UTI, Viral hepatitis |

Abbreviations: Satisfied = Majority criterion satisfied (HR < 1 in ≥ half of basins; directional consistency, not predictive performance); Not satisfied = Majority criterion not satisfied; Asymmetry = mean pairwise transition asymmetry among deceased patients (non-basin states excluded, full observation period); ca. = cancer; conn. = connective; ILD = interstitial lung disease; IHD = ischemic heart disease; Hem. = hematologic; IBD = inflammatory bowel disease; UTI = urinary tract infection; Asp. = aspiration; AF = atrial fibrillation; BPH = benign prostatic hyperplasia; OA = osteoarthritis; RA = rheumatoid arthritis; AKI = acute renal failure; HTN = hypertension; TB = tuberculosis. C numbering is arbitrary and does not imply severity or progression ordering. Asymmetry boundaries (<0.1, 0.1–0.45, ≥0.45) were determined from a gap analysis of the corrected asymmetry distribution (Figure S2). Single-basin diseases have no transition metrics. This classification is descriptive and was not pre-registered or validated on independent data.

```{=openxml}
<w:p><w:pPr><w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/><w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:pPr></w:p>
```

---

## Supplementary figure legends

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S1. Atlas of hospitalization hazard ratios across 51 disease families.

Same axes as Figure 3 (observation density vs. death rate). Point color indicates whether the per-basin majority criterion is satisfied for disease-specific hospitalization (HR < 1 in ≥ half of basins; blue) or not (red); point size reflects |log(HR)| from the hospitalization Cox model. The criterion is satisfied in 35 of 51 disease families.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S2. Asymmetry score versus mean transition rate across multi-basin diseases (C3–C6).

Scatter plot of mean per-step transition rate (x-axis) against mean pairwise asymmetry score (y-axis) for all multi-basin diseases (non-basin states excluded, full observation period). Points are colored by Table 1 category. Natural gaps in the asymmetry distribution (0.02–0.19 and 0.35–0.52) are visible, supporting the three-level classification of multi-basin diseases satisfying the majority criterion into C3 (<0.1), C4 (0.1–0.45), and C5 (≥0.45). Dashed horizontal lines at 0.1 and 0.45 mark the boundaries.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S3. Within-disease basin-to-basin Jaccard similarity.

Dot plot of pairwise Jaccard similarity between basins within each multi-basin disease, computed from the top enriched comorbidity codes per basin. Higher similarity indicates that basins share comorbidity profiles; lower similarity indicates clinically distinct subgroups. Supplements Figure 4.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S4. Cross-disease basin-discriminating comorbidity codes.

Heatmap of the top 10 comorbidity codes with the highest mean basin spread (difference in prevalence between the most- and least-enriched basins) across 33 multi-basin diseases. Values are spread in percentage points. Codes such as chronic periodontitis and hyperlipidemia discriminate basins across many diseases, suggesting shared comorbidity axes. Supplements Figure 4.

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S5. Sepsis — Kaplan–Meier survival curves for all-cause mortality (C6).

Kaplan–Meier survival curves for all-cause mortality in the sepsis cohort, a representative of C6 (multi-basin, majority criterion not satisfied), stratified by basin membership (Basin 1–3; n = 132,559).

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S6. Sepsis — Kaplan–Meier curves for disease-specific hospitalization (C6).

Kaplan–Meier curves for disease-specific hospitalization in the sepsis cohort, stratified by basin membership (Basin 1–3; n = 369,717, excluding admissions before the representative window).

```{=openxml}
<w:p><w:r><w:br w:type="page"/></w:r></w:p>
```

### Figure S7. Sepsis — state occupancy heatmap before death (C6).

State occupancy heatmap before death in the sepsis cohort. Rows are the 25 most frequent Allstates patterns among deceased patients (n = 61,074), sorted by peak timing (states that predominate earlier in the trajectory appear at top; states concentrated near death appear at bottom). Star (★) marks basin states.
