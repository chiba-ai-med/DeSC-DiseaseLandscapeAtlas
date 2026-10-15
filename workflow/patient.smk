include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
               disease=DISEASES),

rule build_patient_scores:
    input:
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        energy="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/E.tsv",
        scores_csv="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.csv",
        scores_mm="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
        person_map="output/{disease}/" + PREPROC + "/cohort/person_map.tsv",
        tekiyo=TEKIYO,
        hosp=HOSP,
    output:
        tsv="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=50000
    benchmark:
        "benchmarks/patient_scores_{disease}.txt"
    log:
        "logs/patient_scores_{disease}.log"
    shell:
        """
        Rscript src/build_patient_scores.R \
            {wildcards.disease} \
            output/{wildcards.disease}/""" + PREPROC + """ \
            {input.tekiyo} {input.hosp} \
            {output.tsv} \
            >& {log}
        """
