# Disease Family Definitions

## Pilot Disease Families

### 1. Dementia (`dementia`)

**Index ICD-10 codes**: F00, F01, F02, F03, G30

- F00: Dementia in Alzheimer's disease
- F01: Vascular dementia
- F02: Dementia in other diseases classified elsewhere
- F03: Unspecified dementia
- G30: Alzheimer's disease

**Rationale**: Dementia patients accumulate downstream complications
(pneumonia, aspiration, UTI, fractures, delirium, malnutrition) that should
form distinct regions in comorbidity space.

**Expected complication codes for severity landscape**:
- J18 (pneumonia), J69 (aspiration pneumonia)
- N39 (UTI)
- S72 (hip fracture), W19 (fall)
- F05 (delirium)
- E46 (malnutrition)
- I50 (heart failure), N18 (CKD)
- I63/I64 (stroke)

### 2. Diabetes Mellitus (`diabetes`)

**Index ICD-10 codes**: E10, E11, E12, E13, E14

- E10: Type 1 diabetes mellitus
- E11: Type 2 diabetes mellitus
- E12: Malnutrition-related diabetes mellitus
- E13: Other specified diabetes mellitus
- E14: Unspecified diabetes mellitus

**Rationale**: Diabetes has a well-characterized progression: microvascular
(retinopathy, nephropathy, neuropathy) then macrovascular (MI, stroke, PAD)
complications. The landscape should distinguish early vs. complicated diabetes.

**Expected complication codes for severity landscape**:
- H36 (diabetic retinopathy)
- N18 (CKD), N08 (diabetic nephropathy)
- G63 (diabetic neuropathy)
- I21 (MI), I63 (stroke)
- I70 (PAD), E78 (hyperlipidemia)
- I10 (hypertension)
- L97 (diabetic foot ulcer)

### 3. Heart Failure / Cancer (`hf_cancer`)

**Index ICD-10 codes**: I50, C00-C97

- I50: Heart failure
- C00-C97: Malignant neoplasms

**Rationale**: These are high-severity conditions where the comorbidity
landscape should naturally reflect end-stage / treatment-related complications.
Combining HF and cancer in one cohort tests whether the landscape can
distinguish different severe-disease trajectories.

**Expected complication codes for severity landscape**:
- Heart failure pathway: I48 (AF), N18 (CKD), J81 (pulmonary edema), I42
  (cardiomyopathy)
- Cancer pathway: D70 (neutropenia), R64 (cachexia), G62 (neuropathy from
  chemo), Z51 (chemo/radio encounters)
- Shared: J18 (pneumonia), N17 (AKI), K72 (liver failure)

## Adding New Disease Families

1. Add a row to `config/disease_families.tsv`
2. Optionally add a `complications` code list in the TSV
3. Re-run `workflow/cohort.smk` onwards

## Code Matching Rules

- ICD-10 codes are matched by **prefix** (e.g., `E10` matches `E10.0`,
  `E10.1`, `E10.9`, etc.)
- The DeSC data uses ICD-10 2013 edition (Japanese modification)
- Codes are stored without dots (e.g., `E109` not `E10.9`) — the matching
  script handles both formats
