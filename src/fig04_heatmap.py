"""Fig04: heatmap only (no dendrogram), ordered by hierarchical clustering.
Two major blocks (Acute/Frailty, Chronic/Lifestyle) + sub-cluster annotations."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering, fcluster
from scipy.spatial.distance import squareform

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

df = pd.read_csv(f"{base}/summary/jaccard_disease_maxjac.tsv", sep='\t', index_col=0)
diseases = df.index.tolist()
mat = df.values
n = len(diseases)

name_map = {
    "dementia": "Dementia", "parkinson": "Parkinson disease",
    "thyroid_disorder": "Thyroid disorder", "diabetes": "Diabetes",
    "anemia": "Anemia", "sepsis": "Sepsis",
    "heart_failure": "Heart failure", "ckd": "CKD",
    "liver_cancer": "Liver cancer", "epilepsy": "Epilepsy",
    "urinary_tract_infection": "UTI", "delirium": "Delirium",
    "aspiration_pneumonia": "Asp. pneumonia",
    "pneumonia": "Pneumonia", "depression": "Depression",
    "asthma": "Asthma", "osteoporosis": "Osteoporosis",
    "osteoarthritis": "Osteoarthritis", "anxiety": "Anxiety",
    "malignancy": "Malignancy", "hypertension": "Hypertension",
    "dyslipidemia": "Dyslipidemia", "sleep_disorders": "Sleep disorders",
    "atrial_fibrillation": "Atrial fibrillation",
    "prostate_cancer": "Prostate cancer",
    "benign_prostatic_hyperplasia": "BPH",
    "rheumatoid_arthritis": "Rheumatoid arthritis",
    "viral_hepatitis": "Viral hepatitis",
    "malnutrition": "Malnutrition",
    "acute_renal_failure": "Acute renal failure",
    "inflammatory_bowel_disease": "IBD",
    "hematologic_malignancy": "Hem. malignancy",
    "peptic_ulcer": "Peptic ulcer",
    "tuberculosis": "Tuberculosis",
    "copd": "COPD", "gout": "Gout", "stroke": "Stroke",
    "ischemic_heart_disease": "IHD",
    "lung_cancer": "Lung cancer",
    "liver_disease": "Liver disease",
    "colorectal_cancer": "Colorectal cancer",
    "stomach_cancer": "Stomach cancer",
    "metastatic_cancer": "Metastatic cancer",
    "pancreatic_cancer": "Pancreatic cancer",
    "breast_cancer": "Breast cancer",
    "obesity": "Obesity",
    "hip_fracture": "Hip fracture",
    "interstitial_lung_disease": "ILD",
    "autoimmune_connective_tissue": "Autoimmune CT",
    "cirrhosis_fibrosis": "Cirrhosis/fibrosis",
}

# Clustering (for ordering only)
dist_condensed = squareform(1 - mat)
Z = linkage(dist_condensed, method='average')
Z = optimal_leaf_ordering(Z, dist_condensed)
leaf_order = leaves_list(Z)

# Reorder: Acute/Frailty at top-left
acute_frailty = {'acute_renal_failure', 'aspiration_pneumonia', 'sepsis'}
ordered = [diseases[i] for i in leaf_order]
acute_ordered = [d for d in ordered if d in acute_frailty]
chronic_ordered = [d for d in ordered if d not in acute_frailty]
final_order = acute_ordered + chronic_ordered

# Build reordered matrix
idx_map = {d: i for i, d in enumerate(diseases)}
final_idx = [idx_map[d] for d in final_order]
heat_mat = mat[np.ix_(final_idx, final_idx)]
heat_names = [name_map.get(d, d) for d in final_order]

# Sub-cluster definitions (contiguous index ranges in final_order, t=0.25)
# Only clusters with 2+ members
subclusters = [
    # (start_idx, end_idx_inclusive, label, color)
    (5,  6,  'Musculoskeletal',       '#2ca02c'),  # green
    (7,  11, 'Immune /\nHematologic', '#9467bd'),  # purple
    (13, 15, 'Stress-related',        '#ff7f0e'),  # orange
    (16, 24, 'Metabolic / CV',        '#1f77b4'),  # blue
    (26, 30, 'GI / Prostate\ncancer', '#8c564b'),  # brown
    (31, 34, 'Cardiovascular',        '#d62728'),  # red
    (35, 37, 'Pulmonary',             '#17becf'),  # teal
    (38, 41, 'Neurodegenerative',     '#7f7f7f'),  # gray
    (43, 44, 'Geriatric acute',       '#bcbd22'),  # olive
]

# ── Figure ──
fig, ax = plt.subplots(figsize=(28, 24))

cmap = LinearSegmentedColormap.from_list('jac',
    [(0, '#2166ac'), (0.25, '#92c5de'), (0.45, '#f7f7f7'),
     (0.65, '#fddbc7'), (0.8, '#f4a582'), (1.0, '#b2182b')])

im = ax.imshow(heat_mat, cmap=cmap, vmin=0.25, vmax=1.0,
               aspect='equal', interpolation='nearest')

ax.set_xticks(range(n))
ax.set_xticklabels(heat_names, rotation=90, fontsize=18, ha='center')
ax.set_yticks(range(n))
ax.set_yticklabels(heat_names, fontsize=18)

ax.set_title('Max Jaccard similarity across basin code sets (top-30 codes)',
             fontsize=24, fontweight='bold', pad=18)

# Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.55, pad=0.02)
cbar.set_label('Jaccard similarity', fontsize=20)
cbar.ax.tick_params(labelsize=16)

# ── Major block annotations ──
n_acute = len(acute_ordered)
n_chronic = len(chronic_ordered)

rect_a = plt.Rectangle((-0.5, -0.5), n_acute, n_acute,
                        linewidth=3, edgecolor='#c0392b', facecolor='none', zorder=10)
ax.add_patch(rect_a)
ax.text(n_acute / 2 - 0.5, -0.9, 'Acute / Frailty',
        fontsize=22, fontweight='bold', color='#c0392b',
        ha='center', va='bottom', zorder=11)

rect_c = plt.Rectangle((n_acute - 0.5, n_acute - 0.5), n_chronic, n_chronic,
                        linewidth=3, edgecolor='#2980b9', facecolor='none', zorder=10)
ax.add_patch(rect_c)
ax.text(n_acute + 2, n_acute - 0.8, 'Chronic / Lifestyle',
        fontsize=22, fontweight='bold', color='#2980b9',
        ha='left', va='bottom', zorder=11)

# ── Sub-cluster annotations ──
for start, end, label, color in subclusters:
    size = end - start + 1
    rect = plt.Rectangle((start - 0.5, start - 0.5), size, size,
                          linewidth=2.5, edgecolor=color, facecolor='none',
                          linestyle=(0, (4, 3)), zorder=10)
    ax.add_patch(rect)
    mid = (start + end) / 2
    # Place label just outside the block, inside the heatmap area
    ax.text(end + 1.5, mid, label,
            fontsize=22.5, fontweight='bold', color=color,
            ha='left', va='center', zorder=11,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      alpha=0.85, edgecolor='none'))

plt.tight_layout()
out = f"{base}/plot/Figures/Fig04/Fig04_jaccard_clustering.png"
fig.savefig(out, dpi=250, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Saved: {out}")
print(f"Acute/Frailty ({n_acute}): {acute_ordered}")
print(f"Chronic/Lifestyle ({n_chronic}): {n_chronic} diseases")
print(f"Sub-clusters annotated: {len(subclusters)}")
