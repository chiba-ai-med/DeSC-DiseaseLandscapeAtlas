"""Regenerate FigS4: cross-disease basin-discriminating codes heatmap.

For each multi-basin disease, compute the "spread" of each ICD-10 code
across basins (max freq_in - min freq_in). Average spread across diseases
to find the top-10 most discriminating codes. Display as a heatmap.

Output: plot/Figures/Supplementary/FigS4_Fig04_crossdisease_discriminating_codes.png
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

name_map = {
    "acute_renal_failure": "Acute renal failure",
    "anemia": "Anemia",
    "anxiety": "Anxiety",
    "aspiration_pneumonia": "Asp. pneumonia",
    "asthma": "Asthma",
    "atrial_fibrillation": "Atrial fibrillation",
    "benign_prostatic_hyperplasia": "BPH",
    "ckd": "CKD",
    "delirium": "Delirium",
    "dementia": "Dementia",
    "depression": "Depression",
    "diabetes": "Diabetes",
    "dyslipidemia": "Dyslipidemia",
    "epilepsy": "Epilepsy",
    "heart_failure": "Heart failure",
    "hematologic_malignancy": "Hem. malignancy",
    "hypertension": "Hypertension",
    "inflammatory_bowel_disease": "IBD",
    "liver_cancer": "Liver ca.",
    "malignancy": "Malignancy",
    "malnutrition": "Malnutrition",
    "osteoarthritis": "Osteoarthritis",
    "osteoporosis": "Osteoporosis",
    "parkinson": "Parkinson",
    "peptic_ulcer": "Peptic ulcer",
    "pneumonia": "Pneumonia",
    "prostate_cancer": "Prostate ca.",
    "rheumatoid_arthritis": "RA",
    "sepsis": "Sepsis",
    "sleep_disorders": "Sleep disorders",
    "thyroid_disorder": "Thyroid disorder",
    "urinary_tract_infection": "UTI",
    "viral_hepatitis": "Viral hepatitis",
}

# Load ICD-10 names
icd10 = pd.read_csv(f"{base}/data/who_icd10_en.tsv", sep="\t")
code_to_name = dict(zip(icd10["icd10_code"], icd10["english_name"]))
code_to_name.setdefault("F009", "Dementia in Alzheimer disease, unspecified")

# Get multi-basin diseases
atlas = pd.read_csv(f"{base}/summary/atlas_summary.tsv", sep="\t")
multi = atlas[(atlas["n_basins"] > 1) & (atlas["outcome"] == "mortality")]
multi = multi[~multi["disease_family"].isin(["hf_cancer", "schizophrenia"])]
diseases = sorted(multi["disease_family"].tolist())

# Compute spread per code per disease
all_spreads = []
for disease in diseases:
    path = f"{base}/plot/{disease}/basin_code_enrichment.tsv"
    df = pd.read_csv(path, sep="\t", quoting=3)
    df["freq_pct"] = df["freq_in"] * 100
    spread = df.groupby("code")["freq_pct"].agg(["max", "min"])
    spread["spread"] = spread["max"] - spread["min"]
    spread = spread.reset_index()
    spread["disease"] = disease
    all_spreads.append(spread[["disease", "code", "spread"]])

spread_df = pd.concat(all_spreads, ignore_index=True)

# Find top-10 codes by mean spread across diseases
mean_spread = spread_df.groupby("code")["spread"].mean().reset_index()
mean_spread.columns = ["code", "mean_spread"]
mean_spread = mean_spread.sort_values("mean_spread", ascending=False)
top10_codes = mean_spread.head(10)["code"].tolist()

# Build disease x code matrix
matrix = spread_df[spread_df["code"].isin(top10_codes)].pivot_table(
    index="disease", columns="code", values="spread", fill_value=0
)
matrix = matrix[top10_codes]

# Human-readable labels
matrix.index = [name_map.get(d, d) for d in matrix.index]
matrix.columns = [code_to_name.get(c, c) for c in matrix.columns]

# Sort diseases by total spread (descending) to match original layout
matrix["_total"] = matrix.sum(axis=1)
matrix = matrix.sort_values("_total", ascending=False)
matrix = matrix.drop(columns="_total")

# Round to integers for annotation
annot = matrix.round(0).astype(int)

fig, ax = plt.subplots(figsize=(22, 18))

sns.heatmap(
    matrix,
    annot=annot,
    fmt="d",
    cmap="Reds",
    vmin=0,
    vmax=70,
    linewidths=0.5,
    linecolor="white",
    ax=ax,
    annot_kws={"size": 16},
    cbar_kws={"label": "Spread\n(pp)", "shrink": 0.6},
)

ax.set_title(
    "Cross-Disease Basin-Discriminating Codes\n"
    f"Top 10 codes by mean basin spread across {len(matrix)} multi-basin diseases",
    fontsize=28,
    pad=20,
)
ax.set_xlabel("", fontsize=24)
ax.set_ylabel("", fontsize=24)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=18)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=16)
cbar.ax.yaxis.label.set_size(18)

plt.tight_layout()
out_path = f"{base}/plot/Figures/Supplementary/FigS4_Fig04_crossdisease_discriminating_codes.png"
fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved: {out_path}")
print(f"Diseases: {len(matrix)}, Codes: {len(top10_codes)}")
print(f"Top 10 codes: {top10_codes}")
