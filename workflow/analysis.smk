include: "common.smk"

rule all:
    input:
        "summary/basin_profiles.tsv",
        "summary/pc_energy_correlation_all.tsv",
        "summary/pc_specificity.tsv",
        "summary/pc_age_sex_energy_correlations.tsv",
        "summary/all_disease_correlation_matrix.csv",
        "summary/per_disease_key_correlations.tsv",

rule basin_pass_fail_profile:
    input:
        expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
               disease=DISEASES),
    output:
        profiles="summary/basin_profiles.tsv",
        scatter="summary/basin_pass_vs_fail_scatter.png",
        death_rate="summary/basin_pass_fail_death_rate.png",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/basin_pass_fail_profile.log"
    shell:
        """
        Rscript src/basin_pass_fail_profile.R >& {log}
        """

rule pc_energy_analysis:
    input:
        expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
               disease=DISEASES),
    output:
        correlations="summary/pc_energy_correlation_all.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/pc_energy_analysis.log"
    shell:
        """
        Rscript src/pc_energy_analysis.R summary >& {log}
        """

rule pc_specificity:
    input:
        "summary/pc_energy_correlation_all.tsv",
    output:
        specificity="summary/pc_specificity.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/pc_specificity.log"
    shell:
        """
        Rscript src/pc_specificity_analysis.R >& {log}
        """

rule pc_axis_correlations:
    input:
        expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
               disease=DISEASES),
    output:
        correlations="summary/pc_age_sex_energy_correlations.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/pc_axis_correlations.log"
    shell:
        """
        Rscript src/pc_axis_correlations.R >& {log}
        """

rule all_disease_correlations:
    input:
        expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
               disease=DISEASES),
    output:
        matrix_csv="summary/all_disease_correlation_matrix.csv",
        per_disease="summary/per_disease_key_correlations.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/all_disease_correlations.log"
    shell:
        """
        Rscript src/all_disease_correlations.R >& {log}
        """
