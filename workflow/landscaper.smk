include: "common.smk"

rule all:
    input:
        expand("output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
               disease=DISEASES),

rule landscaper:
    input:
        scores="output/{disease}/" + PREPROC + "/pca/" + DIM + "/Scores.mm",
    output:
        allstates="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Allstates.tsv",
        basin="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/Basin.tsv",
        subgraph="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/SubGraph.tsv",
        energy="output/{disease}/" + PREPROC + "/landscaper/" + DIM + "/E.tsv",
    threads: 1
    container:
        CONTAINER_LANDSCAPER
    resources:
        mem_mb=100000
    benchmark:
        "benchmarks/landscaper_{disease}.txt"
    log:
        "logs/landscaper_{disease}.log"
    shell:
        """
        outdir=$(dirname {output.allstates})
        tmpdir=$(mktemp -d "./landscaper_tmp.XXXXXX")
        cd "$tmpdir"
        cp /landscaper .
        cp /Snakefile .
        cp -rf /src .
        ./landscaper -i "../{input.scores}" -o "../$outdir" \
            --input_sparse=TRUE --cores=1 --memgb=100
        cd ..
        rm -rf "$tmpdir"
        """
