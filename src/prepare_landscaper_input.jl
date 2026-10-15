## Prepare Landscaper input from PCA Scores
##
## Landscaper expects a Matrix Market file where non-zero entries represent +1
## and absent entries represent -1 in the Ising model. This script verifies
## the Scores.mm format and copies/transforms it if needed.
##
## Usage:
##   julia src/prepare_landscaper_input.jl <scores_mm> <output_mm>

using MatrixMarket
using SparseArrays

scores_path = ARGS[1]
output_path = ARGS[2]

@info "prepare_landscaper_input" input=scores_path output=output_path

mkpath(dirname(output_path))

S = MatrixMarket.mmread(scores_path)
nrows, ncols = size(S)
nnz_count = nnz(S)

@info "Input scores" nrows=nrows ncols=ncols nnz=nnz_count

# Landscaper expects the same format: sparse binary with +1 as non-zero
# Verify all non-zero values are 1 (or truthy)
if eltype(S) <: Bool
    S_out = sparse(S)
elseif eltype(S) <: Integer
    @assert all(x -> x == 1, nonzeros(S)) "Non-binary values in Scores.mm"
    S_out = S
else
    S_out = sparse(S .!= 0)
end

mmwrite(output_path, S_out)
@info "Done" output=output_path
