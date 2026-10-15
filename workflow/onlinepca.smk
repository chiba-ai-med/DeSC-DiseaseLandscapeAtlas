include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
               disease=DISEASES),

rule run_onlinepca:
    input:
        matrix="output/{disease}/" + PREPROC + "/cohort/feature_matrix.coo.zst",
    output:
        scores_csv="output/{disease}/{preproc}/pca/{dim}/Scores.csv",
        scores_mm="output/{disease}/{preproc}/pca/{dim}/Scores.mm",
        eigenvectors="output/{disease}/{preproc}/pca/{dim}/Eigen_vectors.csv",
    threads: 1
    container:
        CONTAINER_MMJL
    resources:
        mem_mb=1000000
    benchmark:
        "benchmarks/onlinepca_{disease}_{preproc}_{dim}.txt"
    log:
        "logs/onlinepca_{disease}_{preproc}_{dim}.log"
    shell:
        """
        mkdir -p /tmp/julia_depot
        export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
        export JULIA_HISTORY="/dev/null"
        /usr/local/julia/bin/julia src/run_onlinepca.jl \
            {input.matrix} \
            {output.scores_csv} \
            {output.scores_mm} \
            {wildcards.dim} \
            >& {log}
        """
