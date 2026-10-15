## Tag cloud of characteristic comorbidity codes per basin
##
## Reads the binary COO matrix, maps to disease codes and basins,
## computes enrichment, and generates tag cloud PNGs.
##
## Usage:
##   julia src/tagcloud_per_basin.jl <disease> <base_dir> <out_dir>

using CairoMakie
using DelimitedFiles
using Printf

disease  = ARGS[1]
base_dir = ARGS[2]
out_dir  = ARGS[3]

mkpath(out_dir)
println("=== Tag cloud per basin: $disease ===")

# --- Load patient_scores.tsv ---
pt_lines = readlines(joinpath(base_dir, "patient", "patient_scores.tsv"))
header = split(pt_lines[1], '\t')
kojin_col = findfirst(==("kojin_id"), header)
basin_col = findfirst(==("basin_id"), header)

kojin_to_basin = Dict{String,Int}()
for i in 2:length(pt_lines)
    fields = split(pt_lines[i], '\t')
    kojin_to_basin[fields[kojin_col]] = parse(Int, fields[basin_col])
end
println("  $(length(kojin_to_basin)) patients")

# --- Load person_map.tsv (take latest per patient) ---
pmap_lines = readlines(joinpath(base_dir, "cohort", "person_map.tsv"))
# new_row_id \t old_row_id \t kojin_id
kojin_best = Dict{String,Tuple{Int,Int}}()  # kojin -> (new_row_id, old_row_id)
for i in 2:length(pmap_lines)
    fields = split(pmap_lines[i], '\t')
    nrid = parse(Int, fields[1])
    orid = parse(Int, fields[2])
    kid = fields[3]
    if !haskey(kojin_best, kid) || orid > kojin_best[kid][2]
        kojin_best[kid] = (nrid, orid)
    end
end

# Map new_row_id -> basin_id (for latest period only)
rowid_to_basin = Dict{Int,Int}()
for (kid, (nrid, _)) in kojin_best
    if haskey(kojin_to_basin, kid)
        rowid_to_basin[nrid] = kojin_to_basin[kid]
    end
end
println("  $(length(rowid_to_basin)) rows mapped to basins")

# --- Load col_ids.txt ---
col_lines = readlines(joinpath(base_dir, "cohort", "col_ids.txt"))
colid_to_code = Dict{Int,String}()
for ln in col_lines
    parts = split(ln, ',')
    colid_to_code[parse(Int, parts[2])] = String(parts[1])
end
n_codes = length(colid_to_code)
println("  $n_codes disease codes")

# --- Read binary COO ---
println("Reading binary COO ...")
coo_path = joinpath(base_dir, "cohort", "feature_matrix.coo.zst")

using CodecZstd, TranscodingStreams
open(coo_path) do fio
    zio = TranscodingStream(ZstdDecompressor(), fio)

    basins = sort(collect(Set(values(rowid_to_basin))))
    n_basins = length(basins)
    println("  Basins: $basins")

    # Count: overall[code] and per_basin[basin][code]
    global overall_count = Dict{String,Int}()
    global basin_count = Dict{Int,Dict{String,Int}}()
    for b in basins
        basin_count[b] = Dict{String,Int}()
    end
    global n_overall = 0
    global n_per_basin = Dict{Int,Int}()
    for b in basins
        n_per_basin[b] = 0
    end

    # Track which rows we've counted (for n_overall / n_per_basin patient count)
    global seen_rows = Set{Int}()

    nnz = 0
    while !eof(zio)
        col_id = read(zio, Int32)
        row_id = read(zio, Int32)
        nnz += 1

        if !haskey(rowid_to_basin, Int(row_id))
            continue
        end

        code = get(colid_to_code, Int(col_id), nothing)
        if code === nothing
            continue
        end

        b = rowid_to_basin[Int(row_id)]
        overall_count[code] = get(overall_count, code, 0) + 1
        basin_count[b][code] = get(basin_count[b], code, 0) + 1

        if !(Int(row_id) in seen_rows)
            push!(seen_rows, Int(row_id))
            n_overall += 1
            n_per_basin[b] = get(n_per_basin, b, 0) + 1
        end
    end
    close(zio)
    println("  $nnz non-zero entries read, $(length(seen_rows)) patients matched")

    # --- Generate tag clouds ---
    for b in basins
        bc = basin_count[b]
        nb = n_per_basin[b]
        if nb < 30
            println("  Basin $b: too few patients ($nb), skipping")
            continue
        end

        # Compute enrichment
        enrichments = Tuple{String,Float64,Float64}[]  # (code, enrichment, freq_in_basin)
        for (code, cnt) in bc
            freq_basin = cnt / nb
            freq_overall = get(overall_count, code, 1) / n_overall
            enrich = freq_basin / freq_overall
            if cnt >= 10
                push!(enrichments, (code, enrich, freq_basin))
            end
        end

        sort!(enrichments, by = x -> -x[2])
        top = enrichments[1:min(50, length(enrichments))]

        if length(top) < 5
            println("  Basin $b: too few codes, skipping")
            continue
        end

        # Plot tag cloud
        fig = Figure(size = (900, 650))
        ax = Axis(fig[1, 1],
            title = "$disease / Basin $b (n=$(@sprintf("%'d", nb)))",
            titlesize = 16)
        hidedecorations!(ax)
        hidespines!(ax)

        # Layout: scatter text with jittered positions
        rng_x = [0.05 + 0.9 * ((i-1) % 10) / 9 for i in 1:length(top)]
        rng_y = [0.95 - 0.9 * ((i-1) ÷ 10) / (length(top) ÷ 10 + 1) for i in 1:length(top)]

        max_enrich = top[1][2]
        min_enrich = top[end][2]

        for (i, (code, enrich, freq)) in enumerate(top)
            fsize = 8 + 20 * (enrich - min_enrich) / max(max_enrich - min_enrich, 1e-6)
            text!(ax, rng_x[i], rng_y[i],
                  text = code,
                  fontsize = fsize,
                  align = (:center, :center),
                  color = (:black, 0.8))
        end

        xlims!(ax, -0.05, 1.05)
        ylims!(ax, -0.05, 1.05)

        out_path = joinpath(out_dir, "tagcloud_basin$(b).png")
        save(out_path, fig, px_per_unit = 2)
        println("  Written: $out_path")
    end

    # --- Save enrichment table ---
    enrich_path = joinpath(out_dir, "basin_code_enrichment.tsv")
    open(enrich_path, "w") do io
        println(io, "basin_id\tcode\tn_basin_code\tfreq_basin\tn_overall_code\tfreq_overall\tenrichment")
        for b in basins
            nb = n_per_basin[b]
            bc = basin_count[b]
            for (code, cnt) in bc
                freq_b = cnt / nb
                freq_o = get(overall_count, code, 1) / n_overall
                enrich = freq_b / freq_o
                @printf(io, "%d\t%s\t%d\t%.6f\t%d\t%.6f\t%.4f\n",
                        b, code, cnt, freq_b,
                        get(overall_count, code, 0), freq_o, enrich)
            end
        end
    end
    println("  Written: $enrich_path")
end

println("Done.")
