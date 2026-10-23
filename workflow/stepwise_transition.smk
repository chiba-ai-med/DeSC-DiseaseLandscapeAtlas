include: "common.smk"

# Only multi-basin diseases
MULTI_BASIN_DISEASES = []
for d in DISEASES:
    bfile = f"output/{d}/all_comorbid/landscaper/7/Basin.tsv"
    import os
    if os.path.exists(bfile):
        with open(bfile) as f:
            nb = len(f.read().strip().split())
        if nb >= 2:
            MULTI_BASIN_DISEASES.append(d)

rule all:
    input:
        expand("output/{disease}/all_comorbid/cox/7/stepwise_transition.tsv",
               disease=MULTI_BASIN_DISEASES),
        "summary/stepwise_transition_summary.tsv"

rule stepwise_transition:
    input:
        traj="plot/{disease}/trajectories_dead.tsv",
        basin="output/{disease}/all_comorbid/landscaper/7/Basin.tsv"
    output:
        "output/{disease}/all_comorbid/cox/7/stepwise_transition.tsv"
    shell:
        "bash src/stepwise_transition_rate.sh {wildcards.disease} {output}"

rule summarize:
    input:
        expand("output/{disease}/all_comorbid/cox/7/stepwise_transition.tsv",
               disease=MULTI_BASIN_DISEASES)
    output:
        "summary/stepwise_transition_summary.tsv"
    run:
        import pandas as pd
        dfs = []
        for f in input:
            df = pd.read_csv(f, sep='\t')
            dfs.append(df)
        all_df = pd.concat(dfs, ignore_index=True)
        # Calculate off-diagonal rate per disease (basin>0 only)
        basin_only = all_df[(all_df['from_basin'] > 0) & (all_df['to_basin'] > 0)]
        summary = []
        for dis in basin_only['disease'].unique():
            d = basin_only[basin_only['disease'] == dis]
            nb = d['n_basins'].iloc[0]
            total = d['count'].sum()
            same = d[d['from_basin'] == d['to_basin']]['count'].sum()
            moved = total - same
            summary.append({
                'disease': dis,
                'n_basins': nb,
                'total_steps': total,
                'same_basin': same,
                'diff_basin': moved,
                'monthly_move_pct': moved / total * 100 if total > 0 else 0
            })
        sdf = pd.DataFrame(summary).sort_values('monthly_move_pct', ascending=False)
        sdf.to_csv(str(output), sep='\t', index=False)
