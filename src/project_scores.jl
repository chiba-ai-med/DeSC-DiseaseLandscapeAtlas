## Project COO data onto pre-computed eigenvectors to get PCA scores.
##
## Avoids re-computing the covariance matrix. Single pass through the COO.
## Outputs Scores.csv and Scores.mm for dim=2..max_dim.
##
## Usage:
##   julia src/project_scores.jl \
##     <coo_zst> <eigenvectors_csv> <colmeanvec_csv> <outdir> <max_dim>

using DelimitedFiles
using SparseArrays
using MatrixMarket

coo_path    = ARGS[1]
eigvec_path = ARGS[2]
colmean_path = ARGS[3]
outdir      = ARGS[4]
max_dim     = parse(Int, ARGS[5])

@info "project_scores" coo=coo_path max_dim=max_dim

# --- Load eigenvectors (n_features x n_features, columns sorted by eigenvalue desc) ---
@info "Loading eigenvectors ..."
V_full = readdlm(eigvec_path, ',', Float64)
n_features = size(V_full, 1)
V = V_full[:, 1:max_dim]  # n_features x max_dim
@info "Eigenvectors" size=size(V)

# --- Load column means ---
colmean = vec(readdlm(colmean_path, ',', Float64))
@assert length(colmean) == n_features "colmean length mismatch"

# --- Read COO (plain binary, decompress .zst on host before calling) ---
io = open(coo_path, "r")
hdr_bytes = read(io, 8)
nrows = reinterpret(Int32, hdr_bytes[1:4])[1]
ncols = reinterpret(Int32, hdr_bytes[5:8])[1]
@info "COO header" nrows=nrows ncols=ncols
@assert ncols == n_features "ncols ($ncols) != n_features ($n_features)"

# --- Stream COO and compute scores ---
# COO is sorted by row. We accumulate entries for each row,
# then project when the row changes.
scores = zeros(Float64, nrows, max_dim)

current_row = Int32(0)
row_cols = Int32[]  # columns with value 1 for current row

function flush_row!(row::Int32, cols::Vector{Int32}, scores::Matrix{Float64},
                    colmean::Vector{Float64}, V::Matrix{Float64})
    row == 0 && return
    n_feat = size(V, 1)
    # Build centered row vector and project
    # x_centered[j] = (1 - colmean[j]) if j in cols, else (0 - colmean[j])
    # scores[row, :] = x_centered' * V
    # Efficient: scores = sum(V[j,:] for j in cols) - sum(colmean .* V, dims=1)
    # Precompute mean_proj = colmean' * V (done outside)
    for j in cols
        @inbounds for d in 1:size(V, 2)
            scores[row, d] += V[j, d]
        end
    end
end

# Precompute mean projection: will subtract after
mean_proj = colmean' * V  # 1 x max_dim

let n_entries = 0, current_row = current_row
    while !eof(io)
        bytes = read(io, 8)
        length(bytes) < 8 && break
        i = reinterpret(Int32, bytes[1:4])[1]
        j = reinterpret(Int32, bytes[5:8])[1]

        if i != current_row
            flush_row!(current_row, row_cols, scores, colmean, V)
            current_row = i
            empty!(row_cols)
        end
        push!(row_cols, j)

        n_entries += 1
        if n_entries % 50_000_000 == 0
            @info "progress" entries=n_entries current_row=current_row
        end
    end
    flush_row!(current_row, row_cols, scores, colmean, V)
    @info "COO read complete" entries=n_entries
end
close(io)

# Subtract mean projection from all rows
for r in 1:nrows
    @inbounds for d in 1:max_dim
        scores[r, d] -= mean_proj[d]
    end
end

# --- Write output for max_dim only ---
dimdir = joinpath(outdir, string(max_dim))
mkpath(dimdir)

writedlm(joinpath(dimdir, "Scores.csv"), scores, ',')

S = sparse(scores .> 0)
mmwrite(joinpath(dimdir, "Scores.mm"), S)

writedlm(joinpath(dimdir, "Eigen_vectors.csv"), V, ',')

@info "Written" dim=max_dim scores_mm=joinpath(dimdir, "Scores.mm")

@info "Done" outdir=outdir
