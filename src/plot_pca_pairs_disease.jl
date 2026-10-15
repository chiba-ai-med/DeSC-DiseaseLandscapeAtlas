## PC1-7 pair density plots per disease, colored by attributes
##
## Usage:
##   /julia_bin/bin/julia src/plot_pca_pairs_disease.jl <disease> <out_dir>
##
## Reads output/<disease>/all_comorbid/patient/patient_scores.tsv

using CairoMakie
using Printf

disease = ARGS[1]
out_dir = ARGS[2]
mkpath(out_dir)

tsv_path = joinpath("output", disease, "all_comorbid/patient/patient_scores.tsv")
@assert isfile(tsv_path) "not found: $tsv_path"

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
const SAMPLE_MAX = 500  # per pair for categorical scatter

pair_index = Dict{Tuple{Int,Int},Int}()
let cnt = 0
    for row in 1:K, col in 1:K
        if row > col
            cnt += 1
            pair_index[(col, row)] = cnt
        end
    end
end

# Label variable definitions
struct LabelVar
    name::String
    tsv_col::String
    kind::Symbol  # :continuous or :categorical
end

const LABEL_VARS = [
    LabelVar("energy",       "energy",          :continuous),
    LabelVar("age",          "age",             :continuous),
    LabelVar("n_admissions", "n_admissions",    :continuous),
    LabelVar("obs_years",    "obs_years",       :continuous),
    LabelVar("death",        "death",           :categorical),
    LabelVar("sex",          "sex",             :categorical),
    LabelVar("kenshin",      "kenshin_data_ari",:categorical),
    LabelVar("hospitalized", "hospitalized",    :categorical),
    LabelVar("insurer",      "insurer",         :categorical),
    LabelVar("basin",        "basin",           :categorical),
]

# ---- Parse header ----
function find_col_indices(header_line::AbstractString)
    cols = split(header_line, '\t')
    d = Dict{String,Int}()
    for (i, c) in enumerate(cols)
        d[String(c)] = i
    end
    return d
end

# ---- Pass 1: min/max of PC1-7 ----
println("[$disease] pass 1: estimating PC min/max ...")
col_min = fill(typemax(Float32), K)
col_max = fill(typemin(Float32), K)
global total_lines = 0
global col_idx_g = Dict{String,Int}()

open(tsv_path, "r") do io
    header = readline(io)
    global col_idx_g = find_col_indices(header)
    for line in eachline(io)
        xs = split(line, '\t')
        for i in 1:K
            ci = col_idx_g["PC$i"]
            v = parse(Float32, xs[ci])
            @inbounds begin
                if v < col_min[i]; col_min[i] = v; end
                if v > col_max[i]; col_max[i] = v; end
            end
        end
        global total_lines += 1
    end
end
println("[$disease] pass 1 done. $total_lines patients")

# ---- Pass 2: build histograms + label stats ----
println("[$disease] pass 2: filling histograms ...")

diag_hists = [zeros(Int, NB) for _ in 1:K]
pair_hists = [zeros(Int, NB, NB) for _ in 1:N_PAIRS]

# Continuous: bin sum/count
struct ContAccum
    pair_sum::Vector{Matrix{Float64}}
    pair_cnt::Vector{Matrix{Int}}
end

struct CatAccum
    sampled_xy::Vector{Vector{Tuple{Float32,Float32}}}
    sampled_labels::Vector{Vector{String}}
    max_per_pair::Int
end

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
            SAMPLE_MAX
        )
    end
end

global lines_done = 0
open(tsv_path, "r") do io
    readline(io)  # skip header
    for line in eachline(io)
        xs = split(line, '\t')

        pcs = Vector{Float32}(undef, K)
        for i in 1:K
            pcs[i] = parse(Float32, xs[col_idx_g["PC$i"]])
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

            for lv in LABEL_VARS
                ci = get(col_idx_g, lv.tsv_col, 0)
                ci == 0 && continue
                raw = String(xs[ci])
                (raw == "" || raw == "NA") && continue

                if lv.kind == :continuous
                    val = tryparse(Float64, raw)
                    val === nothing && continue
                    acc = cont_accums[lv.name]
                    @inbounds acc.pair_sum[idx][bx, by] += val
                    @inbounds acc.pair_cnt[idx][bx, by] += 1
                else
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
        if lines_done % 1_000_000 == 0
            println("  [$disease] $lines_done / $total_lines")
        end
    end
end
println("[$disease] pass 2 done.")

# ---- Drawing functions ----
function draw_density(outfile::String)
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
    Label(fig[0, :], "$disease: PC1-7 density (n=$total_lines)", fontsize=16)
    save(outfile, fig)
    println("  saved: $outfile")
end

function draw_continuous(lv_name::String, outfile::String)
    acc = cont_accums[lv_name]
    pw, ph = 260, 220
    fig = Figure(size = (pw * K + 100, ph * K))
    local hm_ref
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
            S = acc.pair_sum[idx]
            C = acc.pair_cnt[idx]
            M = zeros(Float32, NB, NB)
            for j in 1:NB, i in 1:NB
                M[i, j] = C[i, j] > 0 ? Float32(S[i, j] / C[i, j]) : NaN32
            end
            _set_pretty_ticks!(ax, col_min[col], col_max[col], col_min[row], col_max[row])
            xlims!(ax, (Float64(col_min[col]), Float64(col_max[col])))
            ylims!(ax, (Float64(col_min[row]), Float64(col_max[row])))
            xe = range(col_min[col], col_max[col], length=NB+1)
            ye = range(col_min[row], col_max[row], length=NB+1)
            hm_ref = heatmap!(ax, collect(xe), collect(ye), M;
                              colormap=:plasma, nan_color=:transparent)
        else
            hidespines!(ax); hidedecorations!(ax)
        end
    end
    Colorbar(fig[K, K + 1], hm_ref; label=lv_name)
    Label(fig[0, :], "$disease: PC1-7 by $lv_name (n=$total_lines)", fontsize=16)
    save(outfile, fig)
    println("  saved: $outfile")
end

function draw_categorical(lv_name::String, outfile::String)
    acc = cat_accums[lv_name]
    all_labels = Set{String}()
    for sv in acc.sampled_labels
        for s in sv; push!(all_labels, s); end
    end
    groups = sort(collect(all_labels))
    if isempty(groups)
        println("  skip $lv_name (no data)")
        return
    end
    group2idx = Dict(g => i for (i, g) in enumerate(groups))

    if length(groups) == 2
        gcols = [RGBAf(1.0, 0.2, 0.2, 0.8), RGBAf(0.2, 1.0, 0.2, 0.8)]
    elseif length(groups) <= 10
        tab = Makie.wong_colors()
        gcols = [RGBAf(tab[mod1(i, length(tab))], 0.7) for i in 1:length(groups)]
    else
        gcols = [RGBAf(i / length(groups), 0.3, 1.0 - i / length(groups), 0.6)
                 for i in 1:length(groups)]
    end

    pw, ph = 260, 220
    lw = 180
    fig = Figure(size = (pw * K + lw, ph * K))
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

    elems = [MarkerElement(color=gcols[i], marker=:rect, markersize=12) for i in 1:length(groups)]
    Legend(fig[1, K + 1], elems, groups; title=lv_name)
    Label(fig[0, :], "$disease: PC1-7 by $lv_name (n=$total_lines)", fontsize=16)
    save(outfile, fig)
    println("  saved: $outfile")
end

# ---- Draw all ----
println("\n[$disease] drawing ...")
draw_density(joinpath(out_dir, "pairs_density.png"))

for lv in LABEL_VARS
    if lv.kind == :continuous
        draw_continuous(lv.name, joinpath(out_dir, "pairs_$(lv.name).png"))
    else
        draw_categorical(lv.name, joinpath(out_dir, "pairs_$(lv.name).png"))
    end
end

println("\n[$disease] all done.")
