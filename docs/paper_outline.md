# Paper Outline — Disease Landscape Atlas
Target: Patterns (Cell Press)
Updated: 2026-10-18

## Title (draft)
Energy Landscape Analysis of 51 Disease Families Reveals Cross-Disease
Comorbidity Phenotypes from 11 Million Patient Claims

## Summary (150 words max)
- ELA applied to 51 disease families using DeSC claims (11M patients, 7000+ ICD-10)
- Unsupervised: ICD-10 → PCA(7d) → binarize → Energy Landscape → Basins
- 43/51 diseases PASS Cox sanity check (low energy = high risk)
- Known subtypes auto-rediscovered (renal anemia, cardiorenal syndrome)
- Novel: BasinNA bidirectionality (31 high-risk vs 16 low-risk across diseases)
- Novel: cross-disease frailty/dental signature in high-mortality basins
- Framing: unsupervised phenotyping tool, not prognosis prediction

## Figures

### Fig 1: Analysis Workflow
- Pipeline schematic: DeSC → cohort → feature matrix → PCA → binarize → ELA → Basin
- Show "×51 diseases" scale
- Status: TO BE DRAWN (Illustrator/BioRender)

### Fig 2: PCA Results
- (A) Diabetes PC1-7 pairs density plot
- (B) Thyroid PC1-7 pairs density plot
- (C) TODO: Cumulative variance explained (51 diseases overlaid)
- Message: "PCA captures patient heterogeneity in 7 dimensions"

### Fig 3: Energy Landscape Visualization
- (A) Diabetes landscape summary (basin/energy structure)
- (B) Diabetes KM curves per basin
- (C) Thyroid landscape summary
- (D) Thyroid KM curves per basin
- Message: "Low-energy basins are clinically meaningful stable states"

### Fig 4: 51-Disease Atlas Overview
- (A) Atlas heatmap (51 diseases × key metrics, PASS/FAIL)
- (B) Atlas 3-variable summary
- (C) TODO: Cox sanity check summary + Basin count distribution
- Message: "43/51 diseases PASS; basin structure varies from 1 to 5"

### Fig 5: Validation — Known Subtype Rediscovery
- (A) Anemia Basin1 tagcloud → CKD/PTH codes = renal anemia
- (B) Anemia KM curves
- (C) Heart failure Basin1 tagcloud
- (D) CKD Basin1 tagcloud → HF×CKD Jaccard=0.33 = cardiorenal syndrome
- (E) TODO: Jaccard similarity panel
- Message: "Method rediscovers textbook subtypes unsupervised"

### Fig 6: Diabetes — 4-Basin 2-Axis Structure
- (A-D) Tagclouds for Basin1-4
- (E) TODO: Demographics summary (N, death%, age, sex per basin)
- (F) TODO: 2-axis concept (metabolic axis: B1/B3 vs frailty axis: B2/B4)
- (G) TODO: Ahlqvist comparison concept (orthogonal axes)
- Message: "ELA separates metabolic vs frailty axes; orthogonal to Ahlqvist"

### Fig 7: Thyroid — 5 Basins, Largest Mortality Gap
- (A-D) Tagclouds for Basin1-4
- (E) KM curves showing 20.4% mortality gap (Basin2 vs Basin4)
- (F) TODO: Demographics summary
- Message: "Largest inter-basin mortality gap across 51 diseases; 5 basins unique"

### Fig 8: Sepsis — Acute Disease Basin Structure
- (A-C) Tagclouds for Basin1-3
- (D) KM curves
- (E) TODO: Demographics panel (Basin3: hosp 93.9%, age 82.2)
- Message: "ELA works for acute diseases too; catheter-related infection auto-detected"

### Fig 9: BasinNA Bidirectionality ★★★
- (A) Forest plot (47 diseases, HR of largest basin vs BasinNA)
- (B) TODO: Type A/B concept diagram
- (C) TODO: Representative comparison table (Type A: CKD vs Type B: breast cancer)
- Message: "Non-basin patients are highest-risk in 31 diseases, lowest-risk in 16"

### Fig 10: Cross-Disease Frailty/Dental Signature ★★
- (A) Diabetes Basin4 tagcloud
- (B) Thyroid Basin2 tagcloud
- (C) Malnutrition Basin1 tagcloud
- (D) TODO: Enrichment heatmap (dental/eye codes × 3+ diseases)
- (E) TODO: Diabetes Basin2 vs Basin4 comparison (same codes, 16pt mortality gap)
- Message: "High-mortality basins share dental/eye frailty codes across diseases"

## Sections

### Introduction
- Claims data as untapped resource for disease phenotyping
- Energy Landscape Analysis: from statistical physics to clinical data
- Gap: single-disease ELA exists; no cross-disease atlas
- This paper: 51 disease families, same pipeline, cross-cutting patterns

### Results
1. Pipeline overview + 43/51 PASS Cox sanity check (Fig 1, 4)
2. PCA + Landscape structure (Fig 2, 3)
3. Validation: known subtypes rediscovered (Fig 5)
4. Diabetes: 2-axis structure, Ahlqvist orthogonality (Fig 6)
5. Thyroid: 5 basins, 20% mortality gap (Fig 7)
6. Sepsis: acute disease validation (Fig 8)
7. BasinNA bidirectionality (Fig 9)
8. Cross-disease frailty signature (Fig 10)

### Discussion
- ELA as unsupervised phenotyping (NOT prognosis prediction)
- Cox sanity check = preprocessing validation, not the main claim
- BasinNA bidirectionality: Type A (complex) vs Type B (sparse)
- Frailty signature: oral frailty literature connection
- Ahlqvist orthogonality: claims vs lab-based clustering
- Limitations: ICD-10 coding reliability, no lab values, observational
- Future: prospective validation, integration with clinical data

### STAR Methods
- DeSC data description
- Index cohort definition
- Feature matrix construction (sparse binary COO)
- OnlinePCA (7 dimensions)
- Binarization (±1)
- Energy Landscape Analysis (Landscaper)
- Cox sanity check
- Basin-level enrichment analysis
- Cross-disease similarity (Jaccard)
- Software and containers

## TODO — Remaining Work

### Figure panels to create
- [ ] Fig 1: Full workflow diagram
- [ ] Fig 2C: Cumulative variance explained (51 diseases)
- [ ] Fig 4C: Cox sanity check summary + basin count distribution
- [ ] Fig 5E: HF×CKD Jaccard similarity visualization
- [ ] Fig 6E-G: Diabetes demographics, 2-axis concept, Ahlqvist comparison
- [ ] Fig 7F: Thyroid demographics summary
- [ ] Fig 8E: Sepsis demographics panel
- [ ] Fig 9B-C: BasinNA concept diagram + Type A/B comparison table
- [ ] Fig 10D-E: Enrichment heatmap + Basin2 vs Basin4 comparison

### Analysis to complete
- [ ] Cumulative variance explained data (need to extract from PCA output)
- [ ] HF×CKD enriched code overlap table (beyond Jaccard)
- [ ] Osteoporosis Basin1 CKD-MBD codes (for Fig 5 extension)
- [ ] Update docs/basin_analysis/06_metabolic_endocrine.md Ahlqvist section

### Documentation
- [ ] Update docs with Cox basin membership results (all 10 memos)
- [ ] Remove schizophrenia from BasinNA bidirectionality examples (age confound)
- [ ] Finalize docs/basin_analysis/00_overview.md with figure cross-references

### Writing
- [ ] Draft Summary (150 words)
- [ ] Draft Highlights (3-4 × 85 chars)
- [ ] Introduction
- [ ] Results
- [ ] Discussion
- [ ] STAR Methods
