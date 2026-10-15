## PC1-7 pair density plots (all diseases combined), colored by attributes
##
## Usage:
##   /julia_bin/bin/julia src/plot_pca_pairs_all.jl <out_dir>
##
## Reads output/{disease}/all_comorbid/patient/patient_scores.tsv for all diseases.
## For each color variable, produces a 7×7 lower-triangle density + bin-mean plot.

using CairoMakie
using Printf

# ---- helpers (from old repo Functions.jl) ----
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
const K = 7   # PC1..PC7

out_dir = ARGS[1]
mkpath(out_dir)

# Find all patient_scores.tsv files
base = "output"
diseases = sort(filter(d -> isfile(joinpath(base, d, "all_comorbid/patient/patient_scores.tsv")),
                        readdir(base)))
println("Found $(length(diseases)) diseases")

tsv_paths = [joinpath(base, d, "all_comorbid/patient/patient_scores.tsv") for d in diseases]

# Column definitions in patient_scores.tsv (header-based)
# kojin_id PC1 PC2 PC3 PC4 PC5 PC6 PC7 pattern_id basin_id basin energy age sex death obs_years insurer kenshin_data_ari hospitalized n_admissions
# We parse header once, then use column indices

function find_col_indices(header_line::AbstractString)
    cols = split(header_line, '\t')
    d = Dict{String,Int}()
    for (i, c) in enumerate(cols)
        d[c] = i
    end
    return d
end

# Label variable definitions: name => (tsv_column, kind, parse_func)
# kind: :continuous (bin-mean heatmap) or :categorical (sampled scatter)
struct LabelVar
    name::String
    tsv_col::String
    kind::Symbol    # :continuous or :categorical
end

const LABEL_VARS = [
    LabelVar("energy",       "energy",          :continuous),
    LabelVar("age",          "age",             :continuous),
    LabelVar("death",        "death",           :categorical),
    LabelVar("sex",          "sex",             :categorical),
    LabelVar("kenshin",      "kenshin_data_ari",:categorical),
    LabelVar("hospitalized", "hospitalized",    :categorical),
    LabelVar("n_admissions", "n_admissions",    :continuous),
    LabelVar("obs_years",    "obs_years",       :continuous),
    LabelVar("insurer",      "insurer",         :categorical),
    LabelVar("basin",        "basin",           :categorical),
]

# ---- pair indexing ----
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

# ---- Pass 1: min/max of PC1-7 ----
println("\n[pass 1] estimating PC min/max ...")
col_min = fill(typemax(Float32), K)
col_max = fill(typemin(Float32), K)
global total_lines = 0

for path in tsv_paths
    open(path, "r") do io
        header = readline(io)  # skip header
        for line in eachline(io)
            xs = split(line, '\t')
            for i in 1:K
                v = parse(Float32, xs[i + 1])  # PC1 is column 2
                @inbounds begin
                    if v < col_min[i]; col_min[i] = v; end
                    if v > col_max[i]; col_max[i] = v; end
                end
            end
            global total_lines += 1
        end
    end
end
println("[pass 1] done. total patients = $total_lines")
for i in 1:K
    println("  PC$i: [$(col_min[i]), $(col_max[i])]")
end

# ---- Pass 2: build histograms and accumulate label statistics ----
println("\n[pass 2] building histograms ...")

# Density histograms
diag_hists = [zeros(Int, NB) for _ in 1:K]
pair_hists = [zeros(Int, NB, NB) for _ in 1:N_PAIRS]

# For continuous labels: sum and count per bin for bin-mean
# For categorical labels: sample points
struct ContAccum
    pair_sum::Vector{Matrix{Float64}}   # N_PAIRS matrices (NB x NB)
    pair_cnt::Vector{Matrix{Int}}
end

struct CatAccum
    sampled_xy::Vector{Vector{Tuple{Float32,Float32}}}  # per pair
    sampled_labels::Vector{Vector{String}}               # per pair
    max_per_pair::Int
end

# Initialize accumulators
cont_accums = Dict{String, ContAccum}()
cat_accums = Dict{String, CatAccum}()

for lv in LABEL_VARS
    if lv.kind == :continuous
        cont_accums[lv.name] = ContAccum(
            [zeros(Float64, NB, NB) for _ in 1:N_PAIRS],
            [zeros(Int, NB, NB) for _ in 1:N_PAIRS]
        )
    else
        cat_accums[lv.name] = CatAccum(
            [Tuple{Float32,Float32}[] for _ in 1:N_PAIRS],
            [String[] for _ in 1:N_PAIRS],
            500  # max samples per pair
        )
    end
end

global lines_done = 0
for path in tsv_paths
    open(path, "r") do io
        header_line = readline(io)
        col_idx = find_col_indices(header_line)

        for line in eachline(io)
            xs = split(line, '\t')

            # Parse PC values
            pcs = Vector{Float32}(undef, K)
            for i in 1:K
                pcs[i] = parse(Float32, xs[col_idx["PC$i"]])
            end

            # Diagonal histograms
            for i in 1:K
                b = _safe_bin(pcs[i], col_min[i], col_max[i], NB)
                @inbounds diag_hists[i][b] += 1
            end

            # Pair histograms
            for row in 1:K, col in 1:K
                (row > col) || continue
                bx = _safe_bin(pcs[col], col_min[col], col_max[col], NB)
                by = _safe_bin(pcs[row], col_min[row], col_max[row], NB)
                idx = pair_index[(col, row)]
                @inbounds pair_hists[idx][bx, by] += 1

                # Continuous label accumulation
                for lv in LABEL_VARS
                    if lv.kind == :continuous
                        ci = get(col_idx, lv.tsv_col, 0)
                        ci == 0 && continue
                        raw = xs[ci]
                        (raw == "" || raw == "NA") && continue
                        val = tryparse(Float64, raw)
                        val === nothing && continue
                        acc = cont_accums[lv.name]
                        @inbounds acc.pair_sum[idx][bx, by] += val
                        @inbounds acc.pair_cnt[idx][bx, by] += 1
                    end
                end

                # Categorical label sampling (reservoir-like)
                for lv in LABEL_VARS
                    if lv.kind == :categorical
                        ci = get(col_idx, lv.tsv_col, 0)
                        ci == 0 && continue
                        raw = xs[ci]
                        (raw == "" || raw == "NA") && continue
                        acc = cat_accums[lv.name]
                        if length(acc.sampled_xy[idx]) < acc.max_per_pair
                            if rand() < (acc.max_per_pair / total_lines)
                                push!(acc.sampled_xy[idx], (pcs[col], pcs[row]))
                                push!(acc.sampled_labels[idx], raw)
                            end
                        end
                    end
                end
            end

            global lines_done += 1
            if lines_done % 5_000_000 == 0
                println("  [pass 2] $lines_done / $total_lines")
            end
        end
    end
end
println("[pass 2] done. $lines_done lines processed")

# ---- Drawing ----
println("\n[drawing] ...")

function draw_density_only(outfile::String)
    plot_w = 260
    plot_h = 220
    fig = Figure(size = (plot_w * K, plot_h * K))

    for row in 1:K, col in 1:K
        ax = Axis(fig[row, col])
        if row == col
            xedges = range(col_min[row], col_max[row], length=NB)
            yvals = Float64.(diag_hists[row])
            ymax = maximum(yvals)
            _set_pretty_ticks!(ax, col_min[row], col_max[row], 0.0f0, Float32(ymax))
            lines!(ax, collect(xedges), yvals)
            xlims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            ax.title = "PC$row"
        elseif row > col
            idx = pair_index[(col, row)]
            H = Float32.(pair_hists[idx])
            _set_pretty_ticks!(ax, col_min[col], col_max[col], col_min[row], col_max[row])
            xlims!(ax, (Float64(col_min[col]), Float64(col_max[col])))
            ylims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            xedges = range(col_min[col], col_max[col], length=NB+1)
            yedges = range(col_min[row], col_max[row], length=NB+1)
            heatmap!(ax, collect(xedges), collect(yedges), H; colormap=:viridis)
        else
            hidespines!(ax); hidedecorations!(ax)
        end
    end

    save(outfile, fig)
    println("  saved: $outfile")
end

function draw_continuous_label(lv_name::String, outfile::String)
    acc = cont_accums[lv_name]
    plot_w = 260
    plot_h = 220
    fig = Figure(size = (plot_w * K + 100, plot_h * K))

    for row in 1:K, col in 1:K
        ax = Axis(fig[row, col])
        if row == col
            xedges = range(col_min[row], col_max[row], length=NB)
            yvals = Float64.(diag_hists[row])
            ymax = maximum(yvals)
            _set_pretty_ticks!(ax, col_min[row], col_max[row], 0.0f0, Float32(ymax))
            lines!(ax, collect(xedges), yvals)
            xlims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            ax.title = "PC$row"
        elseif row > col
            idx = pair_index[(col, row)]
            S = acc.pair_sum[idx]
            C = acc.pair_cnt[idx]
            M = zeros(Float32, NB, NB)
            for j in 1:NB, i in 1:NB
                M[i, j] = C[i, j] > 0 ? Float32(S[i, j] / C[i, j]) : NaN32
            end
            _set_pretty_ticks!(ax, col_min[col], col_max[col], col_min[row], col_max[row])
            xlims!(ax, (Float64(col_min[col]), Float64(col_max[col])))
            ylims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            xedges = range(col_min[col], col_max[col], length=NB+1)
            yedges = range(col_min[row], col_max[row], length=NB+1)
            hm = heatmap!(ax, collect(xedges), collect(yedges), M;
                          colormap=:plasma, nan_color=:transparent)
            if row == K && col == 1
                Colorbar(fig[row, K + 1], hm; label=lv_name)
            end
        else
            hidespines!(ax); hidedecorations!(ax)
        end
    end

    Label(fig[0, :], "PC1-7 pairs colored by $lv_name (all diseases, n=$total_lines)",
          fontsize=16)
    save(outfile, fig)
    println("  saved: $outfile")
end

function draw_categorical_label(lv_name::String, outfile::String)
    acc = cat_accums[lv_name]

    # Collect all groups
    all_labels = Set{String}()
    for sv in acc.sampled_labels
        for s in sv; push!(all_labels, s); end
    end
    groups = sort(collect(all_labels))
    group2idx = Dict(g => i for (i, g) in enumerate(groups))

    # Colors
    if length(groups) == 2
        gcols = [RGBAf(1.0, 0.2, 0.2, 0.8), RGBAf(0.2, 1.0, 0.2, 0.8)]
    elseif length(groups) <= 10
        tab10 = Makie.wong_colors()
        gcols = [RGBAf(tab10[mod1(i, length(tab10))], 0.7) for i in 1:length(groups)]
    else
        gcols = [RGBAf(i / length(groups), 0.3, 1.0 - i / length(groups), 0.6)
                 for i in 1:length(groups)]
    end

    plot_w = 260
    plot_h = 220
    legend_w = 180
    fig = Figure(size = (plot_w * K + legend_w, plot_h * K))

    for row in 1:K, col in 1:K
        ax = Axis(fig[row, col])
        if row == col
            xedges = range(col_min[row], col_max[row], length=NB)
            yvals = Float64.(diag_hists[row])
            ymax = maximum(yvals)
            _set_pretty_ticks!(ax, col_min[row], col_max[row], 0.0f0, Float32(ymax))
            lines!(ax, collect(xedges), yvals)
            xlims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            ax.title = "PC$row"
        elseif row > col
            idx = pair_index[(col, row)]
            # Background density
            H = Float32.(pair_hists[idx])
            _set_pretty_ticks!(ax, col_min[col], col_max[col], col_min[row], col_max[row])
            xlims!(ax, (Float64(col_min[col]), Float64(col_max[col])))
            ylims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            xedges = range(col_min[col], col_max[col], length=NB+1)
            yedges = range(col_min[row], col_max[row], length=NB+1)
            heatmap!(ax, collect(xedges), collect(yedges), H; colormap=:viridis)

            # Scatter overlay
            xy = acc.sampled_xy[idx]
            slabels = acc.sampled_labels[idx]
            if !isempty(xy)
                xs = Float32[p[1] for p in xy]
                ys = Float32[p[2] for p in xy]
                pcols = [gcols[get(group2idx, slabels[i], 1)] for i in 1:length(slabels)]
                scatter!(ax, xs, ys; markersize=4, color=pcols)
            end
        else
            hidespines!(ax); hidedecorations!(ax)
        end
    end

    # Legend
    elems = [MarkerElement(color=gcols[i], marker=:rect, markersize=12)
             for i in 1:length(groups)]
    Legend(fig[1, K + 1], elems, groups; title=lv_name)

    Label(fig[0, :], "PC1-7 pairs by $lv_name (all diseases, n=$total_lines)",
          fontsize=16)
    save(outfile, fig)
    println("  saved: $outfile")
end

# Draw density-only
draw_density_only(joinpath(out_dir, "pairs_all_density.png"))

# Draw continuous labels
for lv in LABEL_VARS
    if lv.kind == :continuous
        draw_continuous_label(lv.name, joinpath(out_dir, "pairs_all_$(lv.name).png"))
    end
end

# Draw categorical labels
for lv in LABEL_VARS
    if lv.kind == :categorical
        draw_categorical_label(lv.name, joinpath(out_dir, "pairs_all_$(lv.name).png"))
    end
end

println("\nAll done.")
