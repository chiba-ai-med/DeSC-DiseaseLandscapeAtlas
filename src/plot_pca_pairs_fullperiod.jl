## PC1-7 pair density plots for FULL observation periods (all rows from Scores.csv)
##
## Unlike plot_pca_pairs_disease.jl which reads patient_scores.tsv (latest period only),
## this reads Scores.csv directly (all observation periods, e.g. 21M rows for dementia).
##
## Usage:
##   /julia_bin/bin/julia src/plot_pca_pairs_fullperiod.jl <disease> <out_dir>

using CairoMakie
using Printf

disease = ARGS[1]
out_dir = ARGS[2]
mkpath(out_dir)

scores_path = joinpath("output", disease, "all_comorbid/pca/7/Scores.csv")
@assert isfile(scores_path) "not found: $scores_path"

# ---- helpers ----
@inline function _safe_bin(v::Float32, vmin::Float32, vmax::Float32, nb::Int)
    width = max(vmax - vmin, eps(Float32))
    t = (v - vmin) / width
    b = Int(floor(t * nb)) + 1
    return ifelse(b < 1, 1, ifelse(b > nb, nb, b))
end

@inline function _set_pretty_ticks!(ax::Axis, xmin::Float32, xmax::Float32,
                                     ymin::Float32, ymax::Float32; nticks::Int=3)
    xrange = xmax - xmin
    if abs(xmax) > 1e4 || abs(xmin) > 1e4 || (xrange != 0 && xrange < 1e-3)
        ax.xtickformat = xs -> [@sprintf("%.1e", x) for x in xs]
    end
    yrange = ymax - ymin
    if abs(ymax) > 1e4 || abs(ymin) > 1e4 || (yrange != 0 && yrange < 1e-3)
        ax.ytickformat = ys -> [@sprintf("%.1e", y) for y in ys]
    end
    ax.xticks = Makie.LinearTicks(nticks)
    ax.yticks = Makie.LinearTicks(nticks)
end

# ---- config ----
const NB = 100
const K = 7
const N_PAIRS = K * (K - 1) ÷ 2

pair_index = Dict{Tuple{Int,Int},Int}()
let cnt = 0
    for row in 1:K, col in 1:K
        if row > col
            cnt += 1
            pair_index[(col, row)] = cnt
        end
    end
end

# ---- Pass 1: min/max ----
println("[$disease] pass 1: estimating PC min/max ...")
col_min = fill(typemax(Float32), K)
col_max = fill(typemin(Float32), K)
global total_lines = 0

open(scores_path, "r") do io
    for line in eachline(io)
        xs = split(line, ',')
        for i in 1:K
            v = parse(Float32, xs[i])
            @inbounds begin
                if v < col_min[i]; col_min[i] = v; end
                if v > col_max[i]; col_max[i] = v; end
            end
        end
        global total_lines += 1
    end
end
println("[$disease] pass 1 done. $total_lines rows (all periods)")

# ---- Pass 2: build histograms ----
println("[$disease] pass 2: filling histograms ...")
diag_hists = [zeros(Int, NB) for _ in 1:K]
pair_hists = [zeros(Int, NB, NB) for _ in 1:N_PAIRS]

global lines_done = 0
open(scores_path, "r") do io
    for line in eachline(io)
        xs = split(line, ',')

        pcs = Vector{Float32}(undef, K)
        for i in 1:K
            pcs[i] = parse(Float32, xs[i])
        end

        for i in 1:K
            b = _safe_bin(pcs[i], col_min[i], col_max[i], NB)
            @inbounds diag_hists[i][b] += 1
        end

        for row in 1:K, col in 1:K
            (row > col) || continue
            bx = _safe_bin(pcs[col], col_min[col], col_max[col], NB)
            by = _safe_bin(pcs[row], col_min[row], col_max[row], NB)
            idx = pair_index[(col, row)]
            @inbounds pair_hists[idx][bx, by] += 1
        end

        global lines_done += 1
        if lines_done % 5_000_000 == 0
            println("  [$disease] $lines_done / $total_lines")
        end
    end
end
println("[$disease] pass 2 done.")

# ---- Draw density ----
println("[$disease] drawing ...")
pw, ph = 260, 220
fig = Figure(size = (pw * K, ph * K))
for row in 1:K, col in 1:K
    ax = Axis(fig[row, col])
    if row == col
        xe = range(col_min[row], col_max[row], length=NB)
        yv = Float64.(diag_hists[row])
        _set_pretty_ticks!(ax, col_min[row], col_max[row], 0.0f0, Float32(maximum(yv)))
        lines!(ax, collect(xe), yv)
        xlims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
        ax.title = "PC$row"
    elseif row > col
        idx = pair_index[(col, row)]
        H = Float32.(pair_hists[idx])
        _set_pretty_ticks!(ax, col_min[col], col_max[col], col_min[row], col_max[row])
        xlims!(ax, (Float64(col_min[col]), Float64(col_max[col])))
        ylims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
        xe = range(col_min[col], col_max[col], length=NB+1)
        ye = range(col_min[row], col_max[row], length=NB+1)
        heatmap!(ax, collect(xe), collect(ye), H; colormap=:viridis)
    else
        hidespines!(ax); hidedecorations!(ax)
    end
end
Label(fig[0, :], "$disease: PC1-7 density — all periods (n=$total_lines)", fontsize=16)
outfile = joinpath(out_dir, "pairs_fullperiod.png")
save(outfile, fig, px_per_unit = 2)
println("  saved: $outfile")
println("[$disease] done.")
