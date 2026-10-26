## Cox Energy Sanity Check (optimized)
##
## Vectorized pattern assignment via bit encoding on sparse Scores.mm.
## No dense matrix conversion. Runs in ~3GB peak for 112M rows.
##
## Usage:
##   Rscript src/cox_sanity_check_simple.R \
##     <disease_family> <dim> <base_dir> <tekiyo_csv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
disease   <- args[1]
dim_val   <- as.integer(args[2])
base_dir  <- args[3]
tekiyo_path <- args[4]
out_dir   <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Cox sanity check: %s dim=%d ===\n", disease, dim_val))

# --- Load Landscaper output ---
ls_dir <- file.path(base_dir, "landscaper", dim_val)

sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))
subgraph <- sbmap[, .(pattern_id = state_id, basin_id)]
energy <- sbmap[, .(pattern_id = state_id, energy)]

basin_tbl <- fread(file.path(ls_dir, "Basin.tsv"), header = FALSE)
n_basins <- nrow(basin_tbl)
cat(sprintf("Basins: %d\n", n_basins))

# --- Read Scores.mm (sparse) and assign patterns via bit encoding ---
pca_dir <- file.path(base_dir, "pca", dim_val)
S <- readMM(file.path(pca_dir, "Scores.mm"))
S_t <- as(S, "TsparseMatrix")
n_rows <- nrow(S_t)
k <- ncol(S_t)
n_nz <- length(S_t@i)
cat(sprintf("Scores.mm: %d rows x %d cols (sparse, %d nnz)\n",
            n_rows, k, n_nz))

# Vectorized bit encoding: pattern_key = sum(2^(j-1)) for non-zero columns
i_idx <- S_t@i + 1L  # 1-based row
j_idx <- S_t@j + 1L  # 1-based col
rm(S, S_t)

cat("Computing pattern keys (vectorized) ...\n")
bit_vals <- 2L^(j_idx - 1L)
pattern_keys <- integer(n_rows)  # initialized to 0
# Use data.table for fast grouped sum
ij_dt <- data.table(i = i_idx, bv = bit_vals)
rm(i_idx, j_idx, bit_vals)
key_dt <- ij_dt[, .(key = sum(bv)), by = i]
rm(ij_dt)
pattern_keys[key_dt$i] <- key_dt$key
rm(key_dt)

# Build pattern_key -> pattern_id lookup from Allstates
allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
n_patterns <- nrow(allstates)
cat(sprintf("Allstates: %d patterns\n", n_patterns))

as_keys <- integer(n_patterns)
for (p in seq_len(n_patterns)) {
  as_keys[p] <- sum(2L^(which(allstates[p, ] > 0) - 1L))
}
key2pid <- setNames(seq_len(n_patterns), as_keys)

cat("Mapping pattern keys to pattern IDs ...\n")
unique_keys <- unique(pattern_keys)
key_map <- key2pid[as.character(unique_keys)]
pid_lookup <- setNames(as.integer(key_map), unique_keys)
row_pid <- pid_lookup[as.character(pattern_keys)]
row_pid[is.na(row_pid)] <- 0L
cat(sprintf("Unmatched rows: %d / %d\n", sum(row_pid == 0L), n_rows))
rm(pattern_keys, unique_keys, key_map, pid_lookup)

# --- Build row-level table (lightweight) ---
dt <- data.table(new_row_id = seq_len(n_rows), pattern_id = row_pid)
rm(row_pid)
dt <- merge(dt, subgraph, by = "pattern_id", all.x = TRUE)
dt <- merge(dt, energy, by = "pattern_id", all.x = TRUE)

# --- Load person mapping (only needed columns) ---
person_map_path <- file.path(base_dir, "cohort", "person_map.tsv")
if (!file.exists(person_map_path)) {
  cat("ERROR: person_map.tsv not found\n")
  result <- data.table(
    disease_family = disease, dim = dim_val, outcome = "mortality",
    HR = NA, HR_lower = NA, HR_upper = NA, pvalue = NA,
    direction_ok = NA, n_patients = NA, n_basins = n_basins,
    largest_basin_frac = NA
  )
  fwrite(result, file.path(out_dir, "energy_sanity_check.tsv"), sep = "\t")
  quit(save = "no")
}

cat("Loading person_map ...\n")
pmap <- fread(person_map_path)
dt <- merge(dt, pmap, by = "new_row_id", all.x = TRUE)
rm(pmap)

# --- Take latest row per patient (highest old_row_id) ---
cat("Deduplicating to one row per patient ...\n")
dt_person <- dt[!is.na(kojin_id), .SD[which.max(old_row_id)], by = kojin_id]
rm(dt)
cat(sprintf("Unique patients: %d\n", nrow(dt_person)))

# --- Load tekiyo (only needed columns) ---
cat("Loading tekiyo ...\n")
tek <- fread(tekiyo_path,
             select = c("kojin_id", "birth_ym", "sex_code",
                         "observable_start_ym", "observable_end_ym",
                         "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
dt_person[, kojin_id := as.character(kojin_id)]
dt_person <- merge(dt_person, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

# --- Compute time and events ---
parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
dt_person[, start_y := parse_ym(observable_start_ym)]
dt_person[, end_y := parse_ym(observable_end_ym)]
dt_person[, time_years := end_y - start_y]
dt_person[, event_death := as.integer(!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1)]
dt_person[, age := start_y - parse_ym(birth_ym)]

dt_person <- dt_person[is.finite(time_years) & time_years > 0 &
                         is.finite(energy) & is.finite(age)]
dt_person[, sex_f := factor(sex_code)]
dt_person[, E_z := scale(energy)[, 1]]

cat(sprintf("Valid for Cox: %d patients, %d deaths (%.1f%%)\n",
            nrow(dt_person), sum(dt_person$event_death),
            100 * mean(dt_person$event_death)))

# --- Basin stats ---
basin_freq <- dt_person[, .N, by = basin_id][order(-N)]
largest_basin_frac <- basin_freq$N[1] / nrow(dt_person)
cat("Basin distribution:\n")
print(basin_freq)

# --- Cox: mortality ~ E_z + age + sex ---
if (sum(dt_person$event_death) >= 10) {
  fit <- coxph(Surv(time_years, event_death) ~ E_z + age + sex_f,
               data = dt_person)
  s <- summary(fit)
  hr <- s$conf.int["E_z", ]
  pv <- coef(s)["E_z", 5]

  result <- data.table(
    disease_family = disease, dim = dim_val, outcome = "mortality",
    HR = hr[1], HR_lower = hr[3], HR_upper = hr[4],
    pvalue = pv, direction_ok = (hr[1] < 1.0),
    n_patients = nrow(dt_person), n_basins = n_basins,
    largest_basin_frac = largest_basin_frac
  )
  cat(sprintf("\nMortality: HR(E_z) = %.4f [%.4f, %.4f], p = %.2e -> %s\n",
              hr[1], hr[3], hr[4], pv,
              ifelse(hr[1] < 1, "PASS", "FAIL")))
} else {
  cat("Too few death events\n")
  result <- data.table(
    disease_family = disease, dim = dim_val, outcome = "mortality",
    HR = NA, HR_lower = NA, HR_upper = NA, pvalue = NA, direction_ok = NA,
    n_patients = nrow(dt_person), n_basins = n_basins,
    largest_basin_frac = largest_basin_frac
  )
}

fwrite(result, file.path(out_dir, "energy_sanity_check.tsv"), sep = "\t")
cat(sprintf("Written: %s\n", file.path(out_dir, "energy_sanity_check.tsv")))
