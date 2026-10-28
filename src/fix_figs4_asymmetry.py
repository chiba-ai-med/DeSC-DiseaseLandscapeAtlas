"""
Recompute transition asymmetry from stepwise_transition.tsv (full period, all windows),
excluding basin 0. This fixes the No-basin contamination in the previous computation.

Input: output/{disease}/all_comorbid/cox/7/stepwise_transition.tsv
Output: summary/transition_asymmetry_summary.tsv (overwrite)
        plot/Figures/asymmetry_corrected.tsv (overwrite)
"""
import pandas as pd
import numpy as np
import glob
import os

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

# Table 1 classification
T1 = {
    "Homogeneous": ["breast_cancer","colorectal_cancer","interstitial_lung_disease",
        "ischemic_heart_disease","liver_disease","lung_cancer","metastatic_cancer",
        "obesity","stomach_cancer","stroke"],
    "Reversed": ["autoimmune_connective_tissue","cirrhosis_fibrosis","copd","gout",
        "hip_fracture","pancreatic_cancer","tuberculosis"],
    "Symmetric": ["ckd","hypertension","malnutrition"],
    "Convergent": ["atrial_fibrillation","dementia","heart_failure","sleep_disorders"],
    "Progressive": ["anemia","aspiration_pneumonia","asthma","diabetes",
        "hematologic_malignancy","inflammatory_bowel_disease","malignancy",
        "osteoarthritis","osteoporosis","peptic_ulcer"],
    "Structured reversed": ["acute_renal_failure","anxiety",
        "benign_prostatic_hyperplasia","delirium","depression","dyslipidemia",
        "epilepsy","liver_cancer","parkinson","pneumonia","prostate_cancer",
        "rheumatoid_arthritis","schizophrenia","sepsis","thyroid_disorder",
        "urinary_tract_infection","viral_hepatitis"],
}

def classify(disease):
    for pat, diseases in T1.items():
        if disease in diseases:
            return pat
    return None

files = sorted(glob.glob(f"{base}/output/*/all_comorbid/cox/7/stepwise_transition.tsv"))
print(f"Found {len(files)} stepwise_transition.tsv files")

rows = []
for f in files:
    df = pd.read_csv(f, sep='\t')
    disease = df['disease'].iloc[0]
    n_basins_col = df['n_basins'].iloc[0]  # correct basin count from the file

    # Exclude basin 0 (No basin)
    df_pos = df[(df['from_basin'] > 0) & (df['to_basin'] > 0)]

    # Total steps and movement (basin 0 excluded)
    total_steps = df_pos['count'].sum()
    same_basin = df_pos[df_pos['from_basin'] == df_pos['to_basin']]['count'].sum()
    diff_basin = df_pos[df_pos['from_basin'] != df_pos['to_basin']]['count'].sum()
    move_pct = diff_basin / total_steps * 100 if total_steps > 0 else 0

    # Off-diagonal transitions only
    off = df_pos[df_pos['from_basin'] != df_pos['to_basin']].copy()
    if len(off) == 0:
        continue

    # Compute pairwise asymmetry
    pairs_seen = set()
    asyms = []
    for _, r in off.iterrows():
        i, j = int(r['from_basin']), int(r['to_basin'])
        pair = (min(i, j), max(i, j))
        if pair in pairs_seen:
            continue
        pairs_seen.add(pair)

        r_ij_rows = off[(off['from_basin'] == pair[0]) & (off['to_basin'] == pair[1])]
        r_ji_rows = off[(off['from_basin'] == pair[1]) & (off['to_basin'] == pair[0])]
        rate_ij = r_ij_rows['rate'].values[0] if len(r_ij_rows) > 0 else 0
        rate_ji = r_ji_rows['rate'].values[0] if len(r_ji_rows) > 0 else 0
        denom = rate_ij + rate_ji
        if denom > 0:
            asym = abs(rate_ij - rate_ji) / denom
        else:
            asym = 0
        asyms.append(asym)

    mean_asym = np.mean(asyms) if asyms else np.nan

    rows.append({
        'disease': disease,
        'n_basins': n_basins_col,
        'total_steps': total_steps,
        'same_basin': same_basin,
        'diff_basin': diff_basin,
        'per_step_move_pct': move_pct,
        'mean_asymmetry': mean_asym,
        'n_pairs': len(asyms),
    })

result = pd.DataFrame(rows)
result['pattern'] = result['disease'].apply(classify)
result = result.sort_values('disease')

# Verification checks
print("\n=== VERIFICATION ===")
for d in ['thyroid_disorder', 'anemia', 'ckd', 'dementia']:
    r = result[result['disease'] == d]
    if len(r) > 0:
        r = r.iloc[0]
        print(f"  {d}: n_basins={r['n_basins']}, asymmetry={r['mean_asymmetry']:.4f}, "
              f"n_pairs={r['n_pairs']}, pattern={r['pattern']}")

# Check Table 1 consistency
print("\n=== TABLE 1 CONSISTENCY CHECK ===")
issues = []
for _, r in result.iterrows():
    d = r['disease']
    pat = r['pattern']
    asym = r['mean_asymmetry']
    if pat == 'Symmetric' and asym >= 0.3:
        issues.append(f"  WARNING: {d} classified as Symmetric but asymmetry={asym:.4f} >= 0.3")
    if pat == 'Convergent' and (asym < 0.3 or asym > 0.6):
        issues.append(f"  WARNING: {d} classified as Convergent but asymmetry={asym:.4f} outside [0.3, 0.6]")
    if pat == 'Progressive' and asym < 0.6:
        issues.append(f"  WARNING: {d} classified as Progressive but asymmetry={asym:.4f} < 0.6")

if issues:
    for i in issues:
        print(i)
else:
    print("  All consistent")

# Save
out1 = f"{base}/summary/transition_asymmetry_summary.tsv"
out2 = f"{base}/plot/Figures/asymmetry_corrected.tsv"
result.to_csv(out1, sep='\t', index=False)
result.to_csv(out2, sep='\t', index=False)
print(f"\nSaved: {out1}")
print(f"Saved: {out2}")
print(f"Total diseases: {len(result)}")
