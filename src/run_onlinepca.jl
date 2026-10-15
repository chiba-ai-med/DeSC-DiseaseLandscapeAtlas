## Run OnlinePCA.jl on a cohort feature matrix
##
## Performs exact out-of-core PCA and outputs scores in both CSV and
## Matrix Market (binary +1 only) formats.
##
## Usage:
##   julia src/run_onlinepca.jl <input_coo_zst> <out_scores_csv> <out_scores_mm> <dim>

using OnlinePCA
using MatrixMarket
using SparseArrays

infile      = ARGS[1]
outdir      = dirname(ARGS[2])
outfile_mm  = ARGS[3]
dims        = parse(Int, ARGS[4])

@info "run_onlinepca" input=infile dim=dims

mkpath(outdir)

out = exact_ooc_pca(
    input=infile,
    scale="raw",
    dim=dims,
    chunksize=300000,
    mode="sparse_bincoo"
)

OnlinePCA.output(outdir, out, 0.1f0)

# Save binarized scores as Matrix Market (sparse, +1 entries only)
S = sparse(out[4] .> 0)
mmwrite(outfile_mm, S)

@info "Done" outdir=outdir scores_mm=outfile_mm
