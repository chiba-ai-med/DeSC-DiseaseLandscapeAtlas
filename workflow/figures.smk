include: "common.smk"

rule all:
    input:
        expand("plot/{disease}/" + PREPROC + "/" + DIM + "/landscape_summary.png",
               disease=DISEASES),
        "plot/atlas_heatmap.png",
        expand("output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check_per_basin.tsv",
               disease=DISEASES),
        expand("plot/{disease}/km_per_basin.png",
               disease=DISEASES),
        expand("plot/{disease}/basin_code_enrichment.tsv",
               disease=DISEASES),
        expand("plot/{disease}/state_code_enrichment.tsv",
               disease=DISEASES),
        expand("plot/{disease}/pairs_fullperiod.png",
               disease=DISEASES),
        expand("plot/{disease}/trajectory_heatmap.png",
               disease=DISEASES),
        expand("plot/{disease}/trajectory_energy_per_basin.png",
               disease=DISEASES),
        expand("plot/{disease}/trajectory_entropy.png",
               disease=DISEASES),
        "summary/applicability/disease_metrics.tsv",
        "summary/applicability/entropy_summary.tsv",
        "summary/applicability/panel_B_km_energy_tertile.png",

rule plot_landscape:
    input:
        person="output/{disease}/{preproc}/patient/patient_scores.tsv",
        cox="output/{disease}/{preproc}/cox/{dim}/energy_sanity_check.tsv",
        basin="output/{disease}/{preproc}/landscaper/{dim}/Basin.tsv",
        energy="output/{disease}/{preproc}/landscaper/{dim}/E.tsv",
    output:
        summary_fig="plot/{disease}/{preproc}/{dim}/landscape_summary.png",
        km_fig="plot/{disease}/{preproc}/{dim}/km_energy.png",
        forest_fig="plot/{disease}/{preproc}/{dim}/forest_energy.png",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=10000
    benchmark:
        "benchmarks/plot_{disease}_{preproc}_{dim}.txt"
    log:
        "logs/plot_{disease}_{preproc}_{dim}.log"
    shell:
        """
        Rscript src/plot_landscape_atlas.R \
            {input.person} \
            {input.cox} \
            {input.basin} \
            {input.energy} \
            $(dirname {output.summary_fig}) \
            {wildcards.disease} \
            {wildcards.preproc} \
            {wildcards.dim} \
            >& {log}
        """

rule plot_atlas_heatmap:
    input:
        summary="summary/atlas_summary.tsv",
    output:
        heatmap="plot/atlas_heatmap.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/plot_atlas_heatmap.log"
    shell:
        """
        Rscript src/plot_landscape_atlas.R \
            --atlas-heatmap \
            {input.summary} \
            {output.heatmap} \
            >& {log}
        """

rule cox_per_basin:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        tekiyo=TEKIYO,
    output:
        tsv="output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check_per_basin.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/cox_per_basin_{disease}.log"
    shell:
        """
        Rscript src/cox_per_basin.R \
            {wildcards.disease} \
            {input.scores} \
            {input.tekiyo} \
            $(dirname {output.tsv}) \
            >& {log}
        """

rule km_per_basin:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        tekiyo=TEKIYO,
    output:
        png="plot/{disease}/km_per_basin.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/km_per_basin_{disease}.log"
    shell:
        """
        Rscript src/km_per_basin.R \
            {wildcards.disease} \
            {input.scores} \
            {input.tekiyo} \
            $(dirname {output.png}) \
            >& {log}
        """

rule tagcloud_per_basin:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        coo="output/{disease}/" + PREPROC + "/cohort/feature_matrix.coo.zst",
        col_ids="output/{disease}/" + PREPROC + "/cohort/col_ids.txt",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        scores_mm="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
        icd10_en="data/who_icd10_en.tsv",
    output:
        enrichment="plot/{disease}/basin_code_enrichment.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=50000
    log:
        "logs/tagcloud_{disease}.log"
    shell:
        """
        Rscript src/tagcloud_per_basin.R \
            {wildcards.disease} \
            output/{wildcards.disease}/{PREPROC} \
            $(dirname {output.enrichment}) \
            {input.icd10_en} \
            >& {log}
        """

rule tagcloud_per_state:
    input:
        coo="output/{disease}/" + PREPROC + "/cohort/feature_matrix.coo.zst",
        col_ids="output/{disease}/" + PREPROC + "/cohort/col_ids.txt",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        scores_mm="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
        icd10_en="data/who_icd10_en.tsv",
    output:
        enrichment="plot/{disease}/state_code_enrichment.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=50000
    log:
        "logs/tagcloud_state_{disease}.log"
    shell:
        """
        Rscript src/tagcloud_per_state.R \
            {wildcards.disease} \
            output/{wildcards.disease}/{PREPROC} \
            $(dirname {output.enrichment}) \
            {input.icd10_en} \
            >& {log}
        """

rule trajectory_to_death:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        person_map="output/{disease}/" + PREPROC + "/cohort/person_map.tsv",
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        energy="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/E.tsv",
        scores_mm="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
        window_lookup=config["data"]["window_lookup_path"],
    output:
        heatmap="plot/{disease}/trajectory_heatmap.png",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=80000
    log:
        "logs/trajectory_{disease}.log"
    shell:
        """
        Rscript src/trajectory_to_death.R \
            {wildcards.disease} \
            output/{wildcards.disease}/{PREPROC} \
            $(dirname {output.heatmap}) \
            {input.window_lookup} \
            >& {log}
        """

rule pca_pairs_fullperiod:
    input:
        scores="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.csv",
    output:
        png="plot/{disease}/pairs_fullperiod.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/pairs_fullperiod_{disease}.log"
    shell:
        """
        export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
        export JULIA_HISTORY=/dev/null
        export JULIA_SCRATCH_TRACK_ACCESS=false
        export JULIA_SCRATCH_SPACE=/tmp/julia_scratch
        export TMPDIR=/tmp
        mkdir -p /tmp/julia_depot /tmp/julia_scratch
        /julia_bin/bin/julia src/plot_pca_pairs_fullperiod.jl \
            {wildcards.disease} \
            $(dirname {output.png}) \
            >& {log}
        """

rule trajectory_per_basin:
    input:
        traj="plot/{disease}/trajectories_dead.tsv",
    output:
        energy="plot/{disease}/trajectory_energy_per_basin.png",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=80000
    log:
        "logs/trajectory_per_basin_{disease}.log"
    shell:
        """
        Rscript src/trajectory_per_basin.R \
            {wildcards.disease} \
            $(dirname {output.energy}) \
            >& {log}
        """

rule trajectory_entropy_sankey:
    input:
        traj="plot/{disease}/trajectories_dead.tsv",
    output:
        entropy="plot/{disease}/trajectory_entropy.png",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=80000
    log:
        "logs/trajectory_entropy_sankey_{disease}.log"
    shell:
        """
        Rscript src/trajectory_entropy_sankey.R \
            {wildcards.disease} \
            $(dirname {output.entropy}) \
            >& {log}
        """

rule applicability_analysis:
    input:
        atlas="summary/atlas_summary.tsv",
        patients=expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
                        disease=DISEASES),
        trajectories=expand("plot/{disease}/trajectories_dead.tsv",
                            disease=DISEASES),
    output:
        metrics="summary/applicability/disease_metrics.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/applicability_analysis.log"
    shell:
        """
        Rscript src/applicability_analysis.R >& {log}
        """

rule figure_representative_comparison:
    input:
        atlas="summary/atlas_summary.tsv",
        metrics="summary/applicability/disease_metrics.tsv",
        trajectories=expand("plot/{disease}/trajectories_dead.tsv",
                            disease=["diabetes", "dementia", "heart_failure", "sepsis"]),
    output:
        panel_a="summary/applicability/panel_A_energy_trajectory.png",
        panel_c="summary/applicability/panel_C_state_heatmap.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/figure_representative_comparison.log"
    shell:
        """
        Rscript src/figure_representative_comparison.R >& {log}
        """

rule figure_km_comparison:
    input:
        patients=expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
                        disease=["diabetes", "dementia", "heart_failure", "sepsis"]),
    output:
        panel_b="summary/applicability/panel_B_km_energy_tertile.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/figure_km_comparison.log"
    shell:
        """
        Rscript src/figure_km_comparison.R >& {log}
        """

rule plot_entropy_summary:
    input:
        atlas="summary/atlas_summary.tsv",
        metrics="summary/applicability/disease_metrics.tsv",
        entropies=expand("plot/{disease}/trajectory_entropy.tsv",
                         disease=DISEASES),
        trajectories=expand("plot/{disease}/trajectories_dead.tsv",
                            disease=DISEASES),
    output:
        summary="summary/applicability/entropy_summary.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/plot_entropy_summary.log"
    shell:
        """
        Rscript src/plot_entropy_summary.R >& {log}
        """
