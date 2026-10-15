include: "common.smk"

rule all:
    input:
        expand("output/{disease}/cohort/index_cohort.tsv", disease=DISEASES)

rule define_index_cohort:
    input:
        coo=config["data"]["coo_path"],
        col_ids=config["data"]["col_ids_path"],
        row_ids=config["data"]["row_ids_path"],
        tekiyo=TEKIYO,
        families=config["data"]["disease_families_tsv"],
    output:
        cohort="output/{disease}/cohort/index_cohort.tsv",
    container:
        CONTAINER_ANALYSIS
    resources:
        mem_mb=100000
    benchmark:
        "benchmarks/cohort_{disease}.txt"
    log:
        "logs/cohort_{disease}.log"
    shell:
        """
        Rscript src/define_index_cohort.R \
            {wildcards.disease} \
            {input.families} \
            {input.tekiyo} \
            {input.coo} \
            {input.col_ids} \
            {input.row_ids} \
            {output.cohort} \
            >& {log}
        """
