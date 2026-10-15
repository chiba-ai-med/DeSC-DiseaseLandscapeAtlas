include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/cohort/feature_matrix.coo.zst",
               disease=DISEASES),

rule build_cohort_matrix:
    input:
        cohort="output/{disease}/cohort/index_cohort.tsv",
        coo=config["data"]["coo_path"],
        col_ids=config["data"]["col_ids_path"],
        row_ids=config["data"]["row_ids_path"],
        families=config["data"]["disease_families_tsv"],
    output:
        matrix="output/{disease}/{preproc}/cohort/feature_matrix.coo.zst",
        col_out="output/{disease}/{preproc}/cohort/col_ids.txt",
        row_out="output/{disease}/{preproc}/cohort/row_ids.txt",
        stats="output/{disease}/{preproc}/cohort/matrix_stats.tsv",
    container:
        CONTAINER_MMJL
    resources:
        mem_mb=1000000
    benchmark:
        "benchmarks/matrix_{disease}_{preproc}.txt"
    log:
        "logs/matrix_{disease}_{preproc}.log"
    shell:
        """
        mkdir -p /tmp/julia_depot
        export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
        export JULIA_HISTORY="/dev/null"
        /usr/local/julia/bin/julia src/build_cohort_matrix.jl \
            {wildcards.disease} \
            {wildcards.preproc} \
            {input.cohort} \
            {input.coo} \
            {input.col_ids} \
            {input.row_ids} \
            {input.families} \
            {output.matrix} \
            {output.col_out} \
            {output.row_out} \
            {output.stats} \
            >& {log}
        """
