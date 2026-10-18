include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/comparison/risk_score_comparison.tsv",
               disease=DISEASES),
        "summary/risk_comparison/all_comparisons.tsv",
        "summary/risk_comparison/energy_correlates_summary.tsv",

rule compute_comorbidity_indices:
    input:
        coo="output/{disease}/" + PREPROC + "/cohort/feature_matrix.coo.zst",
        col_ids="output/{disease}/" + PREPROC + "/cohort/col_ids.txt",
        person_map="output/{disease}/" + PREPROC + "/cohort/person_map.tsv",
    output:
        tsv="output/{disease}/" + PREPROC + "/comparison/comorbidity_indices.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=80000
    log:
        "logs/comorbidity_indices_{disease}.log"
    shell:
        """
        Rscript src/compute_comorbidity_indices.R \
            {wildcards.disease} \
            output/{wildcards.disease}/{PREPROC} \
            {output.tsv} \
            >& {log}
        """

rule compare_risk_scores:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        indices="output/{disease}/" + PREPROC + "/comparison/comorbidity_indices.tsv",
        tekiyo=TEKIYO,
        hosp="output/{disease}/" + PREPROC + "/cohort/person_hospitalization_specific.tsv",
    output:
        tsv="output/{disease}/" + PREPROC + "/comparison/risk_score_comparison.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=30000
    log:
        "logs/compare_risk_{disease}.log"
    shell:
        """
        Rscript src/compare_risk_scores.R \
            {wildcards.disease} \
            {input.scores} \
            {input.indices} \
            {input.tekiyo} \
            {input.hosp} \
            $(dirname {output.tsv}) \
            >& {log}
        """

rule summarize_risk_comparison:
    input:
        comparisons=expand("output/{disease}/" + PREPROC + "/comparison/risk_score_comparison.tsv",
                           disease=DISEASES),
    output:
        summary="summary/risk_comparison/all_comparisons.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/summarize_risk_comparison.log"
    shell:
        """
        Rscript src/summarize_risk_comparison.R >& {log}
        """

rule energy_correlates:
    input:
        scores=expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
                       disease=DISEASES),
        indices=expand("output/{disease}/" + PREPROC + "/comparison/comorbidity_indices.tsv",
                        disease=DISEASES),
    output:
        summary="summary/risk_comparison/energy_correlates_summary.tsv",
        detail="summary/risk_comparison/energy_correlates.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=30000
    log:
        "logs/energy_correlates.log"
    shell:
        """
        Rscript src/energy_correlates.R >& {log}
        """
