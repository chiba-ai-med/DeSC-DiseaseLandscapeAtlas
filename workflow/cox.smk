include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check.tsv",
               disease=DISEASES),
        expand("output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check_hosp_specific.tsv",
               disease=DISEASES),

rule cox_energy_sanity_check:
    input:
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        energy="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/E.tsv",
        basin="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Basin.tsv",
        scores="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
        person_map="output/{disease}/" + PREPROC + "/cohort/person_map.tsv",
        tekiyo=TEKIYO,
    output:
        result="output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=50000
    benchmark:
        "benchmarks/cox_{disease}.txt"
    log:
        "logs/cox_{disease}.log"
    shell:
        """
        Rscript src/cox_sanity_check_simple.R \
            {wildcards.disease} """ + DIM + """ \
            output/{wildcards.disease}/""" + PREPROC + """ \
            {input.tekiyo} \
            $(dirname {output.result}) \
            >& {log}
        """

rule cox_hospitalization_specific:
    input:
        scores="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
        hosp="output/{disease}/" + PREPROC + "/cohort/person_hospitalization_specific.tsv",
        tekiyo=TEKIYO,
    output:
        result="output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check_hosp_specific.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=50000
    log:
        "logs/cox_hosp_specific_{disease}.log"
    shell:
        """
        Rscript src/cox_hospitalization_specific.R \
            {wildcards.disease} \
            {input.scores} \
            {input.hosp} \
            {input.tekiyo} \
            {output.result} \
            >& {log}
        """
