## Build patient-level PCA scores + labels table
##
## For each disease family, deduplicates to one row per patient,
## joins with demographics, hospitalization, and landscape assignment.
##
## Usage:
##   Rscript src/build_patient_scores.R <disease> <base_dir> <tekiyo> <hosp> <out_tsv>

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
disease   <- args[1]
base_dir  <- args[2]
tekiyo_path <- args[3]
hosp_path   <- args[4]
out_path    <- args[5]

cat(sprintf("=== build_patient_scores: %s ===\n", disease))

# --- Read Scores.csv (dense, dim=7) ---
scores_path <- file.path(base_dir, "pca", "7", "Scores.csv")
cat("Reading Scores.csv ...\n")
scores <- fread(scores_path, header = FALSE)
setnames(scores, paste0("PC", seq_len(ncol(scores))))
scores[, new_row_id := .I]
n_rows <- nrow(scores)
cat(sprintf("  %d rows x %d PCs\n", n_rows, ncol(scores) - 1))

# --- Load Landscaper assignment via state_basin_map ---
ls_dir <- file.path(base_dir, "landscaper", "7")
map_path <- file.path(ls_dir, "state_basin_map.tsv")
sbmap <- fread(map_path)
subgraph <- sbmap[, .(pattern_id = state_id, basin_id)]
energy_dt <- sbmap[, .(pattern_id = state_id, energy)]

# Assign pattern to each row via bit encoding on Scores.mm
S <- readMM(file.path(base_dir, "pca", "7", "Scores.mm"))
S_t <- as(S, "TsparseMatrix")
i_idx <- S_t@i + 1L
j_idx <- S_t@j + 1L
rm(S, S_t)

ij_dt <- data.table(i = i_idx, bv = 2L^(j_idx - 1L))
rm(i_idx, j_idx)
key_dt <- ij_dt[, .(key = sum(bv)), by = i]
rm(ij_dt)

allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
as_keys <- integer(nrow(allstates))
for (p in seq_len(nrow(allstates))) {
  as_keys[p] <- sum(2L^(which(allstates[p, ] > 0) - 1L))
}
key2pid <- setNames(seq_len(nrow(allstates)), as_keys)

pattern_keys <- integer(n_rows)
pattern_keys[key_dt$i] <- key_dt$key
rm(key_dt)
unique_keys <- unique(pattern_keys)
pid_lookup <- setNames(as.integer(key2pid[as.character(unique_keys)]), unique_keys)
scores[, pattern_id := pid_lookup[as.character(pattern_keys)]]
rm(pattern_keys, pid_lookup)

scores <- merge(scores, subgraph, by = "pattern_id", all.x = TRUE)
scores <- merge(scores, energy_dt, by = "pattern_id", all.x = TRUE)

# --- Load person_map ---
cat("Loading person_map ...\n")
pmap <- fread(file.path(base_dir, "cohort", "person_map.tsv"))
scores <- merge(scores, pmap, by = "new_row_id", all.x = TRUE)
rm(pmap)

# --- Deduplicate: one row per patient (latest period) ---
cat("Deduplicating ...\n")
pt <- scores[!is.na(kojin_id), .SD[which.max(old_row_id)], by = kojin_id]
rm(scores)
cat(sprintf("  %d unique patients\n", nrow(pt)))

# --- Join tekiyo ---
cat("Loading tekiyo ...\n")
tek <- fread(tekiyo_path,
             select = c("kojin_id", "birth_ym", "sex_code",
                         "observable_start_ym", "observable_end_ym",
                         "shibou_flg", "shibou_flg_riyouka",
                         "insurer_shubetsu",
                         "kenshin_data_ari"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

# Compute age, observation period
parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
pt[, age := parse_ym(observable_start_ym) - parse_ym(birth_ym)]
pt[, obs_years := parse_ym(observable_end_ym) - parse_ym(observable_start_ym)]
pt[, death := factor(ifelse(!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1, "dead", "alive"))]
pt[, sex := factor(ifelse(sex_code == 1, "male", "female"))]
pt[, insurer := factor(insurer_shubetsu)]
pt[, basin := factor(paste0("Basin", basin_id))]

# --- Join hospitalization ---
cat("Loading hospitalization ...\n")
hosp <- fread(hosp_path)
hosp[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, hosp, by = "kojin_id", all.x = TRUE)
pt[is.na(n_admissions), n_admissions := 0L]
pt[, hospitalized := factor(ifelse(n_admissions > 0, "yes", "no"))]

# --- Output ---
out_cols <- c("kojin_id", paste0("PC", 1:7),
              "pattern_id", "basin_id", "basin", "energy",
              "age", "sex", "death", "obs_years", "insurer",
              "kenshin_data_ari", "hospitalized", "n_admissions")
out_cols <- intersect(out_cols, names(pt))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(pt[, ..out_cols], out_path, sep = "\t")
cat(sprintf("Written: %s (%d rows, %d cols)\n",
            out_path, nrow(pt), length(out_cols)))
