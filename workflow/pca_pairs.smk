include: "common.smk"

PAIR_VARS = ["density", "energy", "age", "n_admissions", "obs_years",
             "death", "sex", "kenshin", "hospitalized", "insurer", "basin"]

rule all:
    input:
        expand("plot/{disease}/pairs_{var}.png",
               disease=DISEASES, var=PAIR_VARS),

rule pca_pairs:
    input:
        tsv="output/{disease}/" + PREPROC + "/patient/patient_scores.tsv",
    output:
        expand("plot/{{disease}}/pairs_{var}.png", var=PAIR_VARS),
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=16000
    threads: 1
    benchmark:
        "benchmarks/pca_pairs_{disease}.txt"
    log:
        "logs/pca_pairs_{disease}.log"
    shell:
        """
        export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
        export JULIA_HISTORY=/dev/null
        export JULIA_SCRATCH_TRACK_ACCESS=false
        export JULIA_SCRATCH_SPACE=/tmp/julia_scratch
        mkdir -p /tmp/julia_depot /tmp/julia_scratch
        /julia_bin/bin/julia src/plot_pca_pairs_disease.jl \
            {wildcards.disease} plot/{wildcards.disease} \
            >& {log}
        """
