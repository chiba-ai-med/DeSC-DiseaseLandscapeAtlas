"""Regenerate FigS7: sepsis state occupancy heatmap before death, with doubled font sizes."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

base = "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
traj = pd.read_csv(f"{base}/plot/sepsis/trajectories_dead.tsv", sep="\t")

max_dist = 19
hm = traj[traj['periods_before_death'] <= max_dist].copy()

n_at_dist = hm.groupby('periods_before_death')['kojin_id'].nunique().reset_index()
n_at_dist.columns = ['periods_before_death', 'n_patients']

state_dist = hm.groupby(['periods_before_death', 'pattern_id', 'state_label']).size().reset_index(name='n')
state_dist = state_dist.merge(n_at_dist, on='periods_before_death')
state_dist['frac'] = state_dist['n'] / state_dist['n_patients']

state_total = hm.groupby(['pattern_id', 'state_label']).size().reset_index(name='total')
state_total = state_total.sort_values('total', ascending=False)
n_show = min(25, len(state_total))
top_states = state_total.head(n_show).copy()

basin_states = hm[hm['basin_id'] > 0][['pattern_id', 'basin_id']].drop_duplicates()
top_states = top_states.merge(basin_states, on='pattern_id', how='left')
top_states['basin_id'] = top_states['basin_id'].fillna(0).astype(int)
top_states['display_label'] = top_states.apply(
    lambda r: f"★ {r['state_label']}" if r['basin_id'] > 0 else r['state_label'], axis=1)

# Sort: basin states first (by total desc), then non-basin (by total desc)
top_states['sort_key'] = top_states.apply(
    lambda r: (0 if r['basin_id'] > 0 else 1, -r['total']), axis=1)
top_states = top_states.sort_values('sort_key')
label_order = top_states['display_label'].tolist()

plot_data = state_dist[state_dist['pattern_id'].isin(top_states['pattern_id'])]
plot_data = plot_data.merge(
    top_states[['pattern_id', 'display_label']], on='pattern_id')

pivot = plot_data.pivot_table(
    index='display_label', columns='periods_before_death', values='frac', fill_value=0)
pivot = pivot.reindex(label_order)
cols_sorted = sorted(pivot.columns, reverse=True)
pivot = pivot[cols_sorted]

n_dead = hm['kojin_id'].nunique()

fig, ax = plt.subplots(figsize=(16, max(8, n_show * 0.45)))
im = ax.imshow(pivot.values, aspect='auto', cmap='inferno',
               norm=Normalize(vmin=0, vmax=pivot.values.max()))

ax.set_yticks(range(len(label_order)))
ax.set_yticklabels(label_order, fontsize=14)
for tick, label in zip(ax.get_yticklabels(), label_order):
    if label.startswith('★'):
        tick.set_fontweight('bold')
ax.set_xticks(range(len(cols_sorted)))
ax.set_xticklabels([str(c) for c in cols_sorted], fontsize=14)

ax.set_xlabel('Periods before death', fontsize=18)
ax.set_title(f'Sepsis: state occupancy before death (n={n_dead} dead)', fontsize=20)

cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('Fraction\nof patients', fontsize=16)
cbar.ax.tick_params(labelsize=14)

plt.tight_layout()
out = f"{base}/plot/Figures/Supplementary/FigS7_sepsis_trajectory.png"
fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"Saved: {out}")
print(f"  States shown: {n_show}, Dead patients: {n_dead}")
