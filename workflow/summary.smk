include: "common.smk"

rule all:
    input:
        "summary/atlas_summary.tsv"

rule compile_summary:
    input:
        expand("output/{disease}/" + PREPROC + "/cox/" + DIM + "/energy_sanity_check.tsv",
               disease=DISEASES),
    output:
        summary="summary/atlas_summary.tsv",
    container:
        CONTAINER_ANALYSIS
    log:
        "logs/compile_summary.log"
    shell:
        """
        Rscript src/compile_summary.R {output.summary} >& {log}
        """
