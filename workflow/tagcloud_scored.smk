include: "common.smk"

METHODS = ["enrichment", "tfidf", "pmi", "logitdiff", "raw_freq"]

rule all:
    input:
        expand("plot/{disease}/tagcloud_scored/{method}/basin_code_{method}.tsv",
               disease=DISEASES, method=METHODS),

rule tagcloud_scored:
    input:
        enrichment="plot/{disease}/basin_code_enrichment.tsv",
        patient_scores="output/{disease}/all_comorbid/patient/patient_scores.tsv",
    output:
        tsv="plot/{disease}/tagcloud_scored/{method}/basin_code_{method}.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/tagcloud_scored_{disease}_{method}.log"
    shell:
        """
        Rscript src/tagcloud_scored.R \
            {wildcards.disease} \
            {wildcards.method} \
            {input.enrichment} \
            $(dirname {output.tsv}) \
            {input.patient_scores} \
            >& {log}
        """
