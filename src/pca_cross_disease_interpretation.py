#!/usr/bin/env python3
"""
Cross-disease PCA axis interpretation.

Reads saved eigenvectors (Eigen_vectors.csv) and code IDs (col_ids.txt)
for all 50 diseases (excluding hf_cancer and schizophrenia).
For each PC (1-7), extracts top-10 positive and top-10 negative loading
codes per disease, then counts cross-disease frequency.

No PCA recomputation — reads only existing output files.

Output (saved to summary/):
  pc_top_loadings_all.tsv           — per-disease per-PC |loading| top-10
  pc_top_loadings_directional.tsv   — per-disease per-PC direction-specific top-10
  pc_crossdisease_stability.tsv     — |loading|-based, sign consistency
  pc_crossdisease_directional.tsv   — direction-specific cross-disease frequency
  pc_interpretation_table.tsv       — condensed table for paper (5+ diseases)
  pc_best_age_sex_per_disease.tsv   — best age/sex correlated PC per disease

Usage:
  python3 src/pca_cross_disease_interpretation.py
"""
import pandas as pd
import numpy as np
import os
import sys

base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

diseases_df = pd.read_csv(f"{base}/config/disease_families.tsv", sep='\t')
all_diseases = [d for d in diseases_df['family_id'].tolist()
                if d not in ('hf_cancer', 'schizophrenia')]

icd10 = pd.read_csv(f"{base}/data/who_icd10_en.tsv", sep='\t')
code_to_name = dict(zip(icd10['icd10_code'], icd10['english_name']))
code_to_name['M5456'] = 'Low back pain, lumbar region'
code_to_name['M8199'] = 'Osteoporosis, unspecified, site unspecified'
code_to_name['F009'] = 'Dementia in Alzheimer disease, unspecified'

TOP_N = 10
N_PC = 7

# 1. |loading| top-10
rows_abs = []
for disease in all_diseases:
    ev = pd.read_csv(f"{base}/output/{disease}/all_comorbid/pca/7/Eigen_vectors.csv", header=None)
    cols = pd.read_csv(f"{base}/output/{disease}/all_comorbid/cohort/col_ids.txt",
                       header=None, names=['code', 'col_id'])
    codes = cols['code'].values
    for pc_idx in range(N_PC):
        loadings = ev.iloc[:, pc_idx].values
        top_abs = np.argsort(np.abs(loadings))[-TOP_N:][::-1]
        for rank, idx in enumerate(top_abs, 1):
            rows_abs.append({
                'disease': disease, 'pc': pc_idx + 1, 'rank': rank,
                'code': codes[idx], 'loading': loadings[idx],
                'abs_loading': abs(loadings[idx]),
                'sign': '+' if loadings[idx] > 0 else '-',
                'english_name': code_to_name.get(codes[idx], ''),
            })
detail_abs = pd.DataFrame(rows_abs)
detail_abs.to_csv(f"{base}/summary/pc_top_loadings_all.tsv", sep='\t', index=False)

# 2. Direction-specific top-10
rows_dir = []
for disease in all_diseases:
    ev = pd.read_csv(f"{base}/output/{disease}/all_comorbid/pca/7/Eigen_vectors.csv", header=None)
    cols = pd.read_csv(f"{base}/output/{disease}/all_comorbid/cohort/col_ids.txt",
                       header=None, names=['code', 'col_id'])
    codes = cols['code'].values
    for pc_idx in range(N_PC):
        loadings = ev.iloc[:, pc_idx].values
        for direction, indices in [('positive', np.argsort(loadings)[-TOP_N:][::-1]),
                                   ('negative', np.argsort(loadings)[:TOP_N])]:
            for rank, idx in enumerate(indices, 1):
                rows_dir.append({
                    'disease': disease, 'pc': pc_idx + 1,
                    'direction': direction, 'rank': rank,
                    'code': codes[idx], 'loading': loadings[idx],
                    'english_name': code_to_name.get(codes[idx], ''),
                })
detail_dir = pd.DataFrame(rows_dir)
detail_dir.to_csv(f"{base}/summary/pc_top_loadings_directional.tsv", sep='\t', index=False)

# 3. Cross-disease stability (|loading|)
cohort_codes = {}
for d in all_diseases:
    ci = pd.read_csv(f"{base}/output/{d}/all_comorbid/cohort/col_ids.txt",
                     header=None, names=['code', 'col_id'])
    cohort_codes[d] = set(ci['code'].values)

stability_rows = []
for pc in range(1, N_PC + 1):
    pc_data = detail_abs[detail_abs['pc'] == pc]
    for code in pc_data['code'].unique():
        code_data = pc_data[pc_data['code'] == code]
        n_diseases = len(code_data)
        n_pos = (code_data['sign'] == '+').sum()
        n_neg = (code_data['sign'] == '-').sum()
        sign_consistency = max(n_pos, n_neg) / n_diseases if n_diseases > 0 else 0
        n_present = sum(1 for d in all_diseases if code in cohort_codes[d])
        stability_rows.append({
            'pc': pc, 'code': code,
            'english_name': code_to_name.get(code, ''),
            'n_diseases_top10': n_diseases, 'n_diseases_present': n_present,
            'n_positive': n_pos, 'n_negative': n_neg,
            'sign_consistency': round(sign_consistency, 3),
            'dominant_sign': '+' if n_pos >= n_neg else '-',
            'stability': 'stable' if sign_consistency >= 0.75 else 'mixed' if sign_consistency >= 0.6 else 'flip',
        })
stability = pd.DataFrame(stability_rows).sort_values(['pc', 'n_diseases_top10'], ascending=[True, False])
stability.to_csv(f"{base}/summary/pc_crossdisease_stability.tsv", sep='\t', index=False)

# 4. Cross-disease directional summary
dir_summary_rows = []
for pc in range(1, N_PC + 1):
    for direction in ['positive', 'negative']:
        pc_sub = detail_dir[(detail_dir['pc'] == pc) & (detail_dir['direction'] == direction)]
        code_counts = pc_sub.groupby('code')['disease'].nunique().reset_index()
        code_counts.columns = ['code', 'n_diseases']
        for _, row in code_counts.iterrows():
            code = row['code']
            n_present = sum(1 for d in all_diseases if code in cohort_codes[d])
            dir_summary_rows.append({
                'pc': pc, 'direction': direction, 'code': code,
                'english_name': code_to_name.get(code, ''),
                'n_diseases_top10': row['n_diseases'], 'n_diseases_present': n_present,
                'fraction_of_present': round(row['n_diseases'] / max(n_present, 1), 3),
            })
dir_summary = pd.DataFrame(dir_summary_rows).sort_values(['pc', 'direction', 'n_diseases_top10'], ascending=[True, True, False])
dir_summary.to_csv(f"{base}/summary/pc_crossdisease_directional.tsv", sep='\t', index=False)

# 5. Interpretation table (5+ diseases)
interp_rows = []
for pc in range(1, N_PC + 1):
    for direction in ['positive', 'negative']:
        top = dir_summary[(dir_summary['pc'] == pc) & (dir_summary['direction'] == direction)]
        top = top.sort_values('n_diseases_top10', ascending=False).head(10)
        for _, r in top.iterrows():
            if r['n_diseases_top10'] >= 5:
                interp_rows.append({
                    'pc': f"PC{pc}", 'direction': direction,
                    'code': r['code'], 'english_name': r['english_name'],
                    'n_top10': r['n_diseases_top10'], 'n_present': r['n_diseases_present'],
                    'fraction': r['fraction_of_present'],
                })
pd.DataFrame(interp_rows).to_csv(f"{base}/summary/pc_interpretation_table.tsv", sep='\t', index=False)

# 6. Best age/sex PC
cor_df = pd.read_csv(f"{base}/summary/pc_age_sex_energy_correlations.tsv", sep='\t')
cor_df = cor_df[cor_df['pc'].str.startswith('PC') & ~cor_df['disease'].isin(['hf_cancer', 'schizophrenia'])]
age_sex_rows = []
for disease in all_diseases:
    sub = cor_df[cor_df['disease'] == disease]
    if len(sub) == 0:
        continue
    best_age = sub.loc[sub['cor_age'].abs().idxmax()]
    best_sex = sub.loc[sub['cor_sex'].abs().idxmax()]
    age_sex_rows.append({
        'disease': disease,
        'best_age_pc': best_age['pc'], 'age_correlation': round(best_age['cor_age'], 4),
        'best_sex_pc': best_sex['pc'], 'sex_correlation': round(best_sex['cor_sex'], 4),
    })
pd.DataFrame(age_sex_rows).to_csv(f"{base}/summary/pc_best_age_sex_per_disease.tsv", sep='\t', index=False)

print(f"All outputs saved to {base}/summary/")
print(f"  pc_top_loadings_all.tsv ({len(detail_abs)} rows)")
print(f"  pc_top_loadings_directional.tsv ({len(detail_dir)} rows)")
print(f"  pc_crossdisease_stability.tsv ({len(stability)} rows)")
print(f"  pc_crossdisease_directional.tsv ({len(dir_summary)} rows)")
print(f"  pc_interpretation_table.tsv ({len(interp_rows)} rows)")
print(f"  pc_best_age_sex_per_disease.tsv ({len(age_sex_rows)} rows)")
