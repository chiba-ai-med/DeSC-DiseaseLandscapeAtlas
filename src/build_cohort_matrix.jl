## Build feature matrix for a disease-family cohort
##
## Reads the full COO matrix and extracts rows (patients) in the index cohort,
## keeping only columns (disease codes) that pass the preprocessing filter.
## Outputs a new COO matrix in zst format.
##
## Usage:
##   julia src/build_cohort_matrix.jl \
##     <disease_family> <preproc_setting> \
##     <cohort_tsv> <coo_zst> <col_ids> <row_ids> <families_tsv> \
##     <out_matrix> <out_col_ids> <out_row_ids> <out_stats>

using CodecZstd
using TranscodingStreams
using Printf
using DelimitedFiles

disease  = ARGS[1]
preproc  = ARGS[2]
cohort_path   = ARGS[3]
coo_path      = ARGS[4]
col_ids_path  = ARGS[5]
row_ids_path  = ARGS[6]
families_path = ARGS[7]
out_matrix    = ARGS[8]
out_col_ids   = ARGS[9]
out_row_ids   = ARGS[10]
out_stats     = ARGS[11]

@info "build_cohort_matrix" disease=disease preproc=preproc

# --- Load disease family definition ---
families_lines = readlines(families_path)
header = split(strip(families_lines[1]), '\t')
fam_idx = findfirst(==("disease_family"), header)
icd_idx = findfirst(==("index_icd10"), header)
comp_idx = findfirst(==("complication_codes"), header)

index_codes = String[]
complication_codes = String[]
for ln in families_lines[2:end]
    fields = split(strip(ln), '\t')
    if fields[fam_idx] == disease
        index_codes = split(fields[icd_idx], ',')
        if comp_idx !== nothing && length(fields) >= comp_idx
            complication_codes = split(fields[comp_idx], ',')
        end
        break
    end
end
@info "Index codes" codes=index_codes

# --- Load column mapping (code -> col_id) ---
col_lines = readlines(col_ids_path)
col_code = String[]
col_id = Int[]
for ln in col_lines
    ln = replace(ln, "﻿" => "")
    commas = findall(',', ln)
    length(commas) >= 2 || continue
    code = strip(ln[1:commas[1]-1])
    id = parse(Int, strip(ln[commas[1]+1:commas[2]-1]))
    push!(col_code, code)
    push!(col_id, id)
end

# --- Determine which columns to exclude (index disease) ---
function matches_any_prefix(code::String, prefixes::Vector{String})
    any(p -> startswith(code, p), prefixes)
end

is_index = [matches_any_prefix(c, index_codes) for c in col_code]
n_index_cols = sum(is_index)
@info "Columns matching index disease (to exclude)" n=n_index_cols

# --- Load cohort row IDs ---
cohort_lines = readlines(cohort_path)
cohort_header = split(strip(cohort_lines[1]), '\t')
row_id_col = findfirst(==("row_id"), cohort_header)
cohort_row_ids = Set{Int}()
for ln in cohort_lines[2:end]
    fields = split(strip(ln), '\t')
    isempty(fields[row_id_col]) && continue
    push!(cohort_row_ids, parse(Int, fields[row_id_col]))
end
@info "Cohort patients" n=length(cohort_row_ids)

# --- Determine column filter based on preproc setting ---
# First pass: count column frequencies within cohort (excluding index cols)
@info "Pass 1: counting column frequencies in cohort ..."
col_freq = Dict{Int,Int}()
open(coo_path) do fio
    zio = TranscodingStream(ZstdDecompressor(), fio)
    n = 0
    for ln in eachline(zio)
        parts = split(strip(ln))
        i = parse(Int, parts[1])
        j = parse(Int, parts[2])
        if i in cohort_row_ids && !is_index[findfirst(==(j), col_id)]
            col_freq[j] = get(col_freq, j, 0) + 1
        end
        n += 1
        if n % 50_000_000 == 0
            @info "pass 1 progress" lines=n
        end
    end
    close(zio)
end

# Sort columns by frequency (descending)
freq_sorted = sort(collect(col_freq), by = x -> -x[2])

# Build set of columns to keep
keep_cols = Set{Int}()
if preproc == "all_comorbid"
    for (c, _) in freq_sorted
        push!(keep_cols, c)
    end
elseif startswith(preproc, "drop_top")
    n_drop = parse(Int, replace(preproc, "drop_top" => ""))
    drop_set = Set{Int}([c for (c, _) in freq_sorted[1:min(n_drop, length(freq_sorted))]])
    for (c, _) in freq_sorted
        if !(c in drop_set)
            push!(keep_cols, c)
        end
    end
elseif preproc == "complications"
    # Keep only complication-related codes
    for (ci, code) in enumerate(col_code)
        if matches_any_prefix(code, complication_codes) && !is_index[ci]
            push!(keep_cols, col_id[ci])
        end
    end
else
    error("Unknown preproc setting: $preproc")
end

@info "Columns after filtering" n_keep=length(keep_cols) preproc=preproc

# --- Remap columns to contiguous IDs ---
keep_sorted = sort(collect(keep_cols))
old2new = Dict{Int,Int}()
for (new_id, old_id) in enumerate(keep_sorted)
    old2new[old_id] = new_id
end

# Remap rows to contiguous IDs
cohort_rows_sorted = sort(collect(cohort_row_ids))
row_old2new = Dict{Int,Int}()
for (new_id, old_id) in enumerate(cohort_rows_sorted)
    row_old2new[old_id] = new_id
end

# --- Pass 2: write filtered COO ---
@info "Pass 2: writing filtered COO ..."
mkpath(dirname(out_matrix))
nnz = 0
open(out_matrix, "w") do oio
    zout = TranscodingStream(ZstdCompressor(), oio)
    open(coo_path) do fio
        zio = TranscodingStream(ZstdDecompressor(), fio)
        n = 0
        for ln in eachline(zio)
            parts = split(strip(ln))
            i = parse(Int, parts[1])
            j = parse(Int, parts[2])
            if i in cohort_row_ids && j in keep_cols
                new_i = row_old2new[i]
                new_j = old2new[j]
                @printf(zout, "%d %d\n", new_i, new_j)
                nnz += 1
            end
            n += 1
            if n % 50_000_000 == 0
                @info "pass 2 progress" lines=n nnz=nnz
            end
        end
        close(zio)
    end
    close(zout)
end

# --- Write column ID mapping ---
open(out_col_ids, "w") do io
    for old_id in keep_sorted
        ci = findfirst(==(old_id), col_id)
        @printf(io, "%s,%d\n", col_code[ci], old2new[old_id])
    end
end

# --- Write row ID mapping ---
open(out_row_ids, "w") do io
    for old_id in cohort_rows_sorted
        @printf(io, "%d,%d\n", old_id, row_old2new[old_id])
    end
end

# --- Write stats ---
open(out_stats, "w") do io
    println(io, "disease_family\t$disease")
    println(io, "preproc\t$preproc")
    println(io, "n_patients\t$(length(cohort_row_ids))")
    println(io, "n_features\t$(length(keep_cols))")
    println(io, "nnz\t$nnz")
    println(io, "n_index_cols_excluded\t$n_index_cols")
end

@info "Done" out_matrix=out_matrix n_patients=length(cohort_row_ids) n_features=length(keep_cols) nnz=nnz
