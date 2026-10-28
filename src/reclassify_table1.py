"""
Reclassify Table 1 with corrected asymmetry values and redraw FigS4.
Data-driven boundaries from gap analysis:
  Symmetric: asymmetry < 0.1
  Convergent: 0.1 <= asymmetry < 0.45
  Progressive: asymmetry >= 0.45
These replace the old <0.3 / 0.3-0.6 / >0.6 boundaries.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

# Load data
majority = pd.read_csv(f"{base}/plot/Figures/per_disease_majority_summary.tsv", sep='\t')
asym = pd.read_csv(f"{base}/summary/transition_asymmetry_summary.tsv", sep='\t')
atlas = pd.read_csv(f"{base}/summary/atlas_summary.tsv", sep='\t')

# Basin counts from atlas (single-basin diseases have no asym entry)
atlas_mort = atlas[atlas['outcome'] == 'mortality'][['disease_family', 'n_basins', 'direction_ok']].rename(
    columns={'disease_family': 'disease'})

# Merge all
merged = atlas_mort.merge(asym[['disease', 'mean_asymmetry', 'per_step_move_pct', 'n_pairs']],
                          on='disease', how='left')
merged = merged[merged['disease'] != 'hf_cancer']

# Classify
def classify(row):
    d = row['disease']
    nb = row['n_basins']
    cox_ok = row['direction_ok']
    asym_val = row['mean_asymmetry']

    if nb == 1:
        return 'Homogeneous' if cox_ok else 'Reversed'

    # Multi-basin: check per-basin Cox majority
    maj = majority[majority['disease'] == d]
    if len(maj) > 0:
        is_expected = maj.iloc[0]['mort_majority']
    else:
        is_expected = cox_ok

    if not is_expected:
        return 'Structured reversed'

    # Expected multi-basin: classify by asymmetry
    if pd.isna(asym_val):
        return 'Convergent'  # fallback
    if asym_val < 0.1:
        return 'Symmetric'
    elif asym_val < 0.45:
        return 'Convergent'
    else:
        return 'Progressive'

merged['pattern'] = merged.apply(classify, axis=1)

# Build Table 1
pattern_order = ['Homogeneous', 'Reversed', 'Symmetric', 'Convergent', 'Progressive', 'Structured reversed']
name_map = {
    "dementia": "Dementia", "parkinson": "Parkinson disease",
    "thyroid_disorder": "Thyroid disorder", "diabetes": "Diabetes",
    "anemia": "Anemia", "sepsis": "Sepsis",
    "heart_failure": "Heart failure", "ckd": "CKD",
    "liver_cancer": "Liver cancer", "epilepsy": "Epilepsy",
    "urinary_tract_infection": "Urinary tract infection", "delirium": "Delirium",
    "aspiration_pneumonia": "Aspiration pneumonia",
    "pneumonia": "Pneumonia", "depression": "Depression",
    "asthma": "Asthma", "osteoporosis": "Osteoporosis",
    "osteoarthritis": "Osteoarthritis", "anxiety": "Anxiety",
    "malignancy": "Malignancy", "hypertension": "Hypertension",
    "dyslipidemia": "Dyslipidemia", "sleep_disorders": "Sleep disorders",
    "atrial_fibrillation": "Atrial fibrillation",
    "prostate_cancer": "Prostate cancer",
    "benign_prostatic_hyperplasia": "Benign prostatic hyperplasia",
    "rheumatoid_arthritis": "Rheumatoid arthritis",
    "viral_hepatitis": "Viral hepatitis",
    "malnutrition": "Malnutrition",
    "acute_renal_failure": "Acute renal failure",
    "inflammatory_bowel_disease": "Inflammatory bowel disease",
    "hematologic_malignancy": "Hematologic malignancy",
    "peptic_ulcer": "Peptic ulcer",
    "tuberculosis": "Tuberculosis",
    "copd": "COPD", "gout": "Gout", "stroke": "Stroke",
    "ischemic_heart_disease": "Ischemic heart disease",
    "lung_cancer": "Lung cancer",
    "liver_disease": "Liver disease",
    "colorectal_cancer": "Colorectal cancer",
    "stomach_cancer": "Stomach cancer",
    "metastatic_cancer": "Metastatic cancer",
    "pancreatic_cancer": "Pancreatic cancer",
    "breast_cancer": "Breast cancer",
    "obesity": "Obesity",
    "hip_fracture": "Hip fracture",
    "interstitial_lung_disease": "Interstitial lung disease",
    "autoimmune_connective_tissue": "Autoimmune connective tissue",
    "cirrhosis_fibrosis": "Cirrhosis/fibrosis",
    "schizophrenia": "Schizophrenia",
}

table1_rows = []
for pat in pattern_order:
    sub = merged[merged['pattern'] == pat].sort_values('disease')
    diseases_list = sub['disease'].tolist()
    n = len(diseases_list)

    basins = "1" if pat in ['Homogeneous', 'Reversed'] else "Multi"
    cox = "Expected" if pat not in ['Reversed', 'Structured reversed'] else "Reversed"

    if pat == 'Symmetric':
        asym_label = "Low (<0.1)"
    elif pat == 'Convergent':
        asym_label = "Moderate (0.1–0.45)"
    elif pat == 'Progressive':
        asym_label = "High (≥0.45)"
    elif pat == 'Structured reversed':
        asym_label = "Various"
    else:
        asym_label = "-"

    disease_names = ", ".join([name_map.get(d, d) for d in sorted(diseases_list)])
    table1_rows.append({
        'Pattern': pat,
        'Basins': basins,
        'Cox_direction': cox,
        'Asymmetry': asym_label,
        'n': n,
        'Diseases': disease_names,
    })

table1 = pd.DataFrame(table1_rows)
table1_path = f"{base}/plot/Figures/Table1_pattern_classification.tsv"
table1.to_csv(table1_path, sep='\t', index=False)
print(f"Saved: {table1_path}")
print(table1.to_string(index=False))

# Show changes from old classification
old_T1 = {
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
old_classify = {}
for pat, ds in old_T1.items():
    for d in ds:
        old_classify[d] = pat

print("\n=== CHANGES FROM OLD CLASSIFICATION ===")
changes = []
for _, r in merged.iterrows():
    d = r['disease']
    new_pat = r['pattern']
    old_pat = old_classify.get(d, 'Unknown')
    if new_pat != old_pat:
        changes.append(f"  {d:35s} {old_pat:25s} → {new_pat}")

if changes:
    for c in changes:
        print(c)
else:
    print("  No changes")

# ============================================================
# Redraw FigS4 pattern space scatter
# ============================================================
multi = merged[(merged['n_basins'] > 1) & merged['mean_asymmetry'].notna()].copy()
multi = multi[multi['disease'] != 'schizophrenia']
multi['label'] = multi['disease'].map(name_map)

# Only Expected multi-basin for the scatter (Symmetric/Convergent/Progressive)
expected = multi[multi['pattern'].isin(['Symmetric', 'Convergent', 'Progressive'])]
reversed_m = multi[multi['pattern'] == 'Structured reversed']

fig, ax = plt.subplots(figsize=(14, 10))

colors = {
    'Symmetric': '#4DAF4A',
    'Convergent': '#377EB8',
    'Progressive': '#E41A1C',
    'Structured reversed': '#FF7F00',
}
markers = {
    'Symmetric': 's',
    'Convergent': '^',
    'Progressive': 'o',
    'Structured reversed': 'D',
}

for pat in ['Symmetric', 'Convergent', 'Progressive', 'Structured reversed']:
    sub = multi[multi['pattern'] == pat]
    ax.scatter(sub['per_step_move_pct'], sub['mean_asymmetry'],
               c=colors[pat], marker=markers[pat], s=120, alpha=0.85,
               label=pat, edgecolors='white', linewidth=0.5, zorder=5)

# Label all points
for _, r in multi.iterrows():
    ax.annotate(r['label'], (r['per_step_move_pct'], r['mean_asymmetry']),
                fontsize=7, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points')

# Boundary lines
ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.axhline(y=0.45, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.text(ax.get_xlim()[1] * 0.95, 0.05, 'Symmetric', ha='right', va='center',
        fontsize=10, color='gray', style='italic')
ax.text(ax.get_xlim()[1] * 0.95, 0.275, 'Convergent', ha='right', va='center',
        fontsize=10, color='gray', style='italic')
ax.text(ax.get_xlim()[1] * 0.95, 0.65, 'Progressive', ha='right', va='center',
        fontsize=10, color='gray', style='italic')

ax.set_xlabel('Per-step transition rate (%)', fontsize=14)
ax.set_ylabel('Mean transition asymmetry', fontsize=14)
ax.set_title('Multi-basin disease classification\n(basin 0 excluded, full observation period)', fontsize=16)
ax.legend(fontsize=11, loc='upper left')
ax.tick_params(labelsize=11)

plt.tight_layout()
out_path = f"{base}/plot/Figures/Supplementary/FigS_pattern_space.png"
fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"\nSaved: {out_path}")
