"""Regenerate FigS2: asymmetry vs transition rate scatter with C3–C6 labels.

Reads confirmed transition_asymmetry_summary.tsv. Does NOT reclassify;
reads the pattern column as-is and maps old names to C3–C6.

Output: plot/Figures/Supplementary/FigS2_Table1_pattern_classification_space.png
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
os.chdir(base)

df = pd.read_csv("summary/transition_asymmetry_summary.tsv", sep="\t")
# Multi-basin only (same filter as original)
multi = df[df['n_basins'] > 1].copy()
# Exclude schizophrenia (same as original)
multi = multi[multi['disease'] != 'schizophrenia'].copy()

# Map old pattern names to C3–C6
pattern_to_c = {
    'Symmetric': 'C3',
    'Convergent': 'C4',
    'Progressive': 'C5',
    'Structured reversed': 'C6',
}
multi['category'] = multi['pattern'].map(pattern_to_c)

# Human-readable disease labels
name_map = {
    'acute_renal_failure': 'Acute renal failure',
    'anemia': 'Anemia',
    'anxiety': 'Anxiety',
    'aspiration_pneumonia': 'Aspiration pneumonia',
    'asthma': 'Asthma',
    'atrial_fibrillation': 'Atrial fibrillation',
    'benign_prostatic_hyperplasia': 'Benign prostatic hyperplasia',
    'ckd': 'CKD',
    'delirium': 'Delirium',
    'dementia': 'Dementia',
    'depression': 'Depression',
    'diabetes': 'Diabetes',
    'dyslipidemia': 'Dyslipidemia',
    'epilepsy': 'Epilepsy',
    'heart_failure': 'Heart failure',
    'hematologic_malignancy': 'Hematologic malignancy',
    'hypertension': 'Hypertension',
    'inflammatory_bowel_disease': 'Inflammatory bowel disease',
    'liver_cancer': 'Liver cancer',
    'malignancy': 'Malignancy',
    'malnutrition': 'Malnutrition',
    'osteoarthritis': 'Osteoarthritis',
    'osteoporosis': 'Osteoporosis',
    'parkinson': 'Parkinson disease',
    'peptic_ulcer': 'Peptic ulcer',
    'pneumonia': 'Pneumonia',
    'prostate_cancer': 'Prostate cancer',
    'rheumatoid_arthritis': 'Rheumatoid arthritis',
    'sepsis': 'Sepsis',
    'sleep_disorders': 'Sleep disorders',
    'thyroid_disorder': 'Thyroid disorder',
    'urinary_tract_infection': 'Urinary tract infection',
    'viral_hepatitis': 'Viral hepatitis',
}
multi['label'] = multi['disease'].map(name_map).fillna(multi['disease'])

# Colors and markers matching original (same hues, new labels)
colors = {
    'C3': '#4DAF4A',   # green (was Symmetric)
    'C4': '#377EB8',   # blue (was Convergent)
    'C5': '#E41A1C',   # red (was Progressive)
    'C6': '#FF7F00',   # orange (was Structured reversed)
}
markers = {
    'C3': 's',   # square
    'C4': '^',   # triangle
    'C5': 'o',   # circle
    'C6': 'D',   # diamond
}

fig, ax = plt.subplots(figsize=(14, 10))

for cat in ['C3', 'C4', 'C5', 'C6']:
    sub = multi[multi['category'] == cat]
    ax.scatter(sub['per_step_move_pct'], sub['mean_asymmetry'],
               c=colors[cat], marker=markers[cat], s=120, alpha=0.85,
               label=cat, edgecolors='white', linewidth=0.5, zorder=5)

# Labels
for _, r in multi.iterrows():
    ax.annotate(r['label'], (r['per_step_move_pct'], r['mean_asymmetry']),
                fontsize=14, ha='left', va='bottom',
                xytext=(5, 5), textcoords='offset points')

# Boundary lines
ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.axhline(y=0.45, color='gray', linestyle='--', alpha=0.5, linewidth=1)

# Region labels on right margin (C6 is not a region — it spans all levels)
ax.text(ax.get_xlim()[1] * 0.95, 0.05, 'C3', ha='right', va='center',
        fontsize=20, color='gray', style='italic')
ax.text(ax.get_xlim()[1] * 0.95, 0.275, 'C4', ha='right', va='center',
        fontsize=20, color='gray', style='italic')
ax.text(ax.get_xlim()[1] * 0.95, 0.65, 'C5', ha='right', va='center',
        fontsize=20, color='gray', style='italic')

ax.set_xlabel('Per-step transition rate (%)', fontsize=22)
ax.set_ylabel('Mean transition asymmetry', fontsize=22)
ax.set_title('Multi-basin disease classification\n(non-basin states excluded, full observation period)', fontsize=24)
ax.legend(fontsize=18, loc='lower left')
ax.tick_params(labelsize=18)

plt.tight_layout()
out_path = "plot/Figures/Supplementary/FigS2_pattern_space.png"
fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
out_pdf = out_path.replace('.png', '.pdf')
fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"Saved: {out_path}")
print(f"Saved: {out_pdf}")

# Verification
print(f"\nVerification:")
print(f"  Diseases plotted: {len(multi)} (expected 33 = 34 multi-basin minus schizophrenia)")
for cat in ['C3', 'C4', 'C5', 'C6']:
    n = len(multi[multi['category'] == cat])
    print(f"  {cat}: {n} diseases")
