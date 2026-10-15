## Re-draw basin pair plots for ALL diseases in a single Julia session
## (avoids repeated CairoMakie JIT compilation)

using CairoMakie
using Printf

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

const NB = 100
const K = 7
const N_PAIRS = K * (K - 1) ÷ 2
const SAMPLE_MAX = 500

pair_index = Dict{Tuple{Int,Int},Int}()
let cnt = 0
    for row in 1:K, col in 1:K
        if row > col
            cnt += 1
            pair_index[(col, row)] = cnt
        end
    end
end

function find_col_indices(header_line::AbstractString)
    cols = split(header_line, '\t')
    d = Dict{String,Int}()
    for (i, c) in enumerate(cols)
        d[String(c)] = i
    end
    return d
end

function process_disease(disease::String)
    tsv_path = joinpath("output", disease, "all_comorbid/patient/patient_scores.tsv")
    out_dir = joinpath("plot", disease)
    isfile(tsv_path) || return println("  SKIP $disease")
    mkpath(out_dir)

    # Pass 1: min/max
    col_min = fill(typemax(Float32), K)
    col_max = fill(typemin(Float32), K)
    total = 0
    local cidx::Dict{String,Int}

    open(tsv_path, "r") do io
        cidx = find_col_indices(readline(io))
        for line in eachline(io)
            xs = split(line, '\t')
            for i in 1:K
                v = parse(Float32, xs[cidx["PC$i"]])
                if v < col_min[i]; col_min[i] = v; end
                if v > col_max[i]; col_max[i] = v; end
            end
            total += 1
        end
    end

    # Pass 2: density + basin sampling
    diag_hists = [zeros(Int, NB) for _ in 1:K]
    pair_hists = [zeros(Int, NB, NB) for _ in 1:N_PAIRS]
    sampled_xy = [Tuple{Float32,Float32}[] for _ in 1:N_PAIRS]
    sampled_labels = [String[] for _ in 1:N_PAIRS]

    open(tsv_path, "r") do io
        readline(io)
        for line in eachline(io)
            xs = split(line, '\t')
            pcs = Vector{Float32}(undef, K)
            for i in 1:K
                pcs[i] = parse(Float32, xs[cidx["PC$i"]])
            end
            for i in 1:K
                b = _safe_bin(pcs[i], col_min[i], col_max[i], NB)
                @inbounds diag_hists[i][b] += 1
            end
            basin_ci = get(cidx, "basin", 0)
            basin_val = basin_ci > 0 ? String(xs[basin_ci]) : ""
            for row in 1:K, col in 1:K
                (row > col) || continue
                bx = _safe_bin(pcs[col], col_min[col], col_max[col], NB)
                by = _safe_bin(pcs[row], col_min[row], col_max[row], NB)
                idx = pair_index[(col, row)]
                @inbounds pair_hists[idx][bx, by] += 1
                if basin_val != "" && basin_val != "NA"
                    if length(sampled_xy[idx]) < SAMPLE_MAX && rand() < (SAMPLE_MAX / total)
                        push!(sampled_xy[idx], (pcs[col], pcs[row]))
                        push!(sampled_labels[idx], basin_val)
                    end
                end
            end
        end
    end

    # Draw
    all_labels = Set{String}()
    for sv in sampled_labels; for s in sv; push!(all_labels, s); end; end
    groups = sort(collect(all_labels))
    isempty(groups) && return println("  $disease: no basin data")
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

    pw, ph, lw = 260, 220, 180
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
            xy = sampled_xy[idx]
            sl = sampled_labels[idx]
            if !isempty(xy)
                xs = Float32[p[1] for p in xy]
                ys = Float32[p[2] for p in xy]
                pc = [gcols[get(group2idx, sl[i], 1)] for i in 1:length(sl)]
                scatter!(ax, xs, ys; markersize=4, color=pc)
            end
        else
            hidespines!(ax); hidedecorations!(ax)
        end
    end
    elems = [MarkerElement(color=gcols[i], marker=:rect, markersize=12) for i in 1:length(groups)]
    Legend(fig[1, K + 1], elems, groups; title="basin")
    Label(fig[0, :], "$disease: PC1-7 by basin (n=$total)", fontsize=16)
    save(joinpath(out_dir, "pairs_basin.png"), fig)
    println("  $disease: pairs_basin.png saved")
end

# Main loop
diseases = sort(filter(d -> isfile(joinpath("output", d, "all_comorbid/patient/patient_scores.tsv")),
                        readdir("output")))
println("Re-drawing basin pair plots for $(length(diseases)) diseases")
for disease in diseases
    process_disease(disease)
end
println("All done.")
