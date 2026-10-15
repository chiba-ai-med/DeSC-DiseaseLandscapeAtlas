# Energy Cox Sanity Check Policy

## Purpose

The Energy Cox sanity check is a **go/no-go gate** for each analysis setting.
It determines whether the resulting landscape has the expected directionality
before investing in downstream interpretation.

## Models

For each `{disease_family} x {preproc} x {dim}`:

```
Model M: Surv(time, death)          ~ E_z + age + sex
Model H: Surv(time, hospitalization) ~ E_z + age + sex
```

where `E_z` is the within-cohort Z-standardized energy from the Landscaper.

## Expected Direction

We expect landscapes where **low-energy basins capture clinically severe /
post-complication states**. In the Cox model, this means:

| Metric              | Expected         | Interpretation                       |
|---------------------|------------------|--------------------------------------|
| HR per +1 SD Energy | < 1.0            | Lower energy -> higher risk          |
| energy_direction_ok | TRUE             | Landscape reflects disease severity  |

## Decision Rules

### Pass (energy_direction_ok = TRUE)

- HR(E_z) < 1.0 for at least one of mortality or hospitalization
- Proceed with full downstream analysis (basin interpretation, transition
  analysis, etc.)

### Fail (energy_direction_ok = FALSE)

- HR(E_z) >= 1.0 for both mortality and hospitalization
- The landscape is a **frequency landscape** (atypical recording patterns =
  high risk)
- Do NOT proceed with downstream interpretation as a severity landscape
- Consider:
  1. Trying a different preprocessing setting
  2. Adding biomarker features
  3. Revisiting the cohort definition

### Borderline

- HR(E_z) close to 1.0 (e.g., 0.95-1.05) with wide CI
- Report but flag as inconclusive
- May need larger cohort or different dimensionality

## Reporting

Each setting produces `output/{disease}/{preproc}/cox/{dim}/energy_sanity_check.tsv`:

```
outcome	HR	HR_lower	HR_upper	pvalue	direction_ok
mortality	0.82	0.75	0.90	0.0001	TRUE
hospitalization	0.88	0.82	0.95	0.001	TRUE
```

The summary table aggregates these across all settings for comparison.

## Important Notes

- This is a **sanity check**, not a final result for publication
- The goal is to screen preprocessing strategies efficiently
- A passing sanity check is necessary but not sufficient for clinical
  interpretability — basin annotation and domain knowledge are still required
- Energy is standardized **within each cohort** (not globally) to enable
  cross-disease comparison of HRs
