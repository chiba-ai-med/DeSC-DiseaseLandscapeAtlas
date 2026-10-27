## Landmark Cox for all 51 diseases using compact binary window lookup
##
## Usage: Rscript landmark_cox_all_diseases.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)
SCRATCH <- "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad"
DATA_START <- "2014/04"
DATA_END   <- "2022/02"

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}

decode_ym <- function(encoded) {
  yr <- encoded %/% 12L
  mo <- encoded %% 12L
  sprintf("%04d/%02d", yr, mo)
}

add_months_str <- function(ym_str, n) {
  y <- as.integer(substr(ym_str, 1, 4))
  m <- as.integer(substr(ym_str, 6, 7))
  m2 <- m + n
  y2 <- y + (m2 - 1L) %/% 12L
  m2 <- ((m2 - 1L) %% 12L) + 1L
  sprintf("%04d/%02d", y2, m2)
}

# ============================================================
# STEP 1: Load compact binary window lookup (1.4GB, all 347M rows)
# ============================================================
cat("=== Loading compact window lookup ===\n")
bin_path <- file.path(SCRATCH, "window_lookup.bin")
stopifnot(file.exists(bin_path))

con <- file(bin_path, "rb")
raw <- readBin(con, "integer", n = 347764079L * 2L, size = 2L, signed = FALSE, endian = "little")
close(con)
cat(sprintf("  Loaded %d uint16 values (%.1f GB)\n", length(raw), length(raw) * 2 / 1e9))

# Reshape: raw[1] = start of row 1, raw[2] = end of row 1, raw[3] = start of row 2, ...
# Index: start[i] = raw[2*(i-1)+1], end[i] = raw[2*(i-1)+2]
# For vectorized access: use even/odd indexing
win_start_enc <- raw[seq(1L, length(raw), 2L)]  # odd positions = starts
win_end_enc   <- raw[seq(2L, length(raw), 2L)]  # even positions = ends
rm(raw)
gc()
cat(sprintf("  Window lookup ready: %d rows\n", length(win_start_enc)))

# ============================================================
# STEP 2: Load tekiyo (once for all diseases)
# ============================================================
cat("\n=== Loading tekiyo ===\n")
tek <- fread("data/tekiyo.csv",
             select = c("kojin_id", "birth_ym", "sex_code",
                        "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
cat(sprintf("  %d patients\n", nrow(tek)))

# ============================================================
# STEP 3: Get disease list
# ============================================================
diseases <- fread("config/disease_families.tsv", sep = "\t")
disease_list <- diseases[family_id != "hf_cancer", family_id]
cat(sprintf("\n=== Processing %d diseases ===\n\n", length(disease_list)))

# ============================================================
# STEP 4: Process each disease
# ============================================================
overall_results <- list()
basin_results <- list()

for (di in seq_along(disease_list)) {
  disease <- disease_list[di]
  base_dir <- file.path("output", disease, "all_comorbid")
  cat(sprintf("[%d/%d] %s ...\n", di, length(disease_list), disease))
  t0 <- proc.time()["elapsed"]

  # --- Load person_map ---
  pmap_path <- file.path(base_dir, "cohort", "person_map.tsv")
  pmap <- fread(pmap_path)
  pmap[, kojin_id := as.character(kojin_id)]
  n_obs <- nrow(pmap)
  n_patients <- uniqueN(pmap$kojin_id)

  # --- Window lookup via binary index ---
  valid_idx <- pmap$old_row_id >= 1L & pmap$old_row_id <= length(win_start_enc)
  pmap[valid_idx, `:=`(
    win_s_enc = win_start_enc[old_row_id],
    win_e_enc = win_end_enc[old_row_id]
  )]
  pmap[valid_idx & win_s_enc > 0L, window_start := decode_ym(win_s_enc)]
  pmap[valid_idx & win_e_enc > 0L, window_end := decode_ym(win_e_enc)]
  pmap[, c("win_s_enc", "win_e_enc") := NULL]

  # --- Landscape assignment ---
  ls_dir <- file.path(base_dir, "landscaper", "7")
  sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))

  scores_mm_path <- file.path(base_dir, "pca", "7", "Scores.mm")
  if (!file.exists(scores_mm_path)) {
    cat(sprintf("  SKIP: no Scores.mm\n"))
    next
  }

  # Read MM header to get nrows, then stream COO entries via fread
  hdr <- readLines(scores_mm_path, n = 2)
  dims <- as.integer(strsplit(trimws(hdr[2]), " +")[[1]])
  n_rows_ls <- dims[1]

  ij_dt <- fread(scores_mm_path, skip = 2, header = FALSE,
                 col.names = c("i", "j"), select = 1:2)
  ij_dt[, bv := 2L^(j - 1L)]
  key_dt <- ij_dt[, .(key = sum(bv)), by = i]
  rm(ij_dt)

  allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
  as_keys <- integer(nrow(allstates))
  for (p in seq_len(nrow(allstates))) {
    as_keys[p] <- sum(2L^(which(allstates[p, ] > 0) - 1L))
  }
  key2pid <- setNames(seq_len(nrow(allstates)), as_keys)

  # Map bit-key → pattern_id directly on key_dt (avoids n_rows_ls-length vector)
  key_dt[, pattern_id := as.integer(key2pid[as.character(key)])]
  setnames(key_dt, "i", "new_row_id")

  # Merge with state_basin_map to get basin_id and energy
  key_dt <- merge(key_dt, sbmap[, .(pattern_id = state_id, basin_id, energy)],
                  by = "pattern_id", all.x = TRUE)

  pmap <- merge(pmap, key_dt[, .(new_row_id, basin_id, energy)],
                by = "new_row_id", all.x = TRUE)
  rm(key_dt)

  # Fill all-negative state rows (absent from Scores.mm → NA after merge)
  all_neg_pid <- as.integer(key2pid["0"])
  if (!is.na(all_neg_pid)) {
    neg_info <- sbmap[state_id == all_neg_pid]
    if (nrow(neg_info) > 0) {
      pmap[is.na(basin_id), `:=`(basin_id = neg_info$basin_id[1],
                                  energy = neg_info$energy[1])]
    }
  }

  # --- Merge tekiyo ---
  pmap <- merge(pmap, tek, by = "kojin_id", all.x = TRUE)

  # --- Eligibility ---
  pmap[, eligible := (
    !is.na(window_start) &
    window_start >= observable_start_ym &
    window_end <= observable_end_ym &
    window_start >= DATA_START &
    window_end <= DATA_END &
    !is.na(basin_id)
  )]

  n_eligible_pts <- uniqueN(pmap[eligible == TRUE, kojin_id])

  # --- First eligible window ---
  eligible_obs <- pmap[eligible == TRUE]
  if (nrow(eligible_obs) == 0) {
    cat(sprintf("  SKIP: no eligible observations\n"))
    overall_results[[disease]] <- data.table(
      disease = disease, design = "landmark",
      HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = 0L, deaths = 0L, direction = "SKIP"
    )
    rm(pmap, eligible_obs); gc()
    next
  }

  setorder(eligible_obs, kojin_id, old_row_id)
  landmark <- eligible_obs[!duplicated(kojin_id)]
  rm(eligible_obs)

  # --- Landmark timing ---
  landmark[, landmark_ym := add_months_str(window_end, 1L)]
  landmark[, landmark_y := parse_ym(landmark_ym)]
  landmark[, obs_start_y := parse_ym(observable_start_ym)]
  landmark[, obs_end_y := parse_ym(observable_end_ym)]
  landmark[, followup_time := obs_end_y - landmark_y]

  # --- Exclusions ---
  landmark[, riyouka := as.integer(shibou_flg_riyouka)]
  landmark[, dead := (!is.na(shibou_flg) & shibou_flg == 1)]
  landmark[, died_before := (dead & obs_end_y <= landmark_y)]
  landmark[, left_before := (!dead & obs_end_y <= landmark_y)]

  cox_pop <- landmark[riyouka == 1 &
                      !died_before &
                      !left_before &
                      followup_time > 0]

  if (nrow(cox_pop) < 100 || sum(cox_pop$dead) < 20) {
    cat(sprintf("  SKIP: too few (N=%d, deaths=%d)\n", nrow(cox_pop), sum(cox_pop$dead)))
    overall_results[[disease]] <- data.table(
      disease = disease, design = "landmark",
      HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = nrow(cox_pop), deaths = sum(cox_pop$dead), direction = "SKIP"
    )
    rm(pmap, landmark, cox_pop); gc()
    next
  }

  cox_pop[, birth_y := parse_ym(birth_ym)]
  cox_pop[, age_at_landmark := landmark_y - birth_y]
  cox_pop[, event := as.integer(dead)]
  cox_pop[, E_z := scale(energy)[,1]]

  n_final <- nrow(cox_pop)
  n_deaths <- sum(cox_pop$event)
  n_excluded_riy <- sum(landmark$riyouka == 0, na.rm = TRUE)
  n_died_before <- sum(landmark$died_before, na.rm = TRUE)
  n_left_before <- sum(landmark$left_before, na.rm = TRUE)

  # --- Overall Cox ---
  fit <- tryCatch(
    coxph(Surv(followup_time, event) ~ E_z + age_at_landmark + factor(sex_code),
          data = cox_pop),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    cat(sprintf("  MODEL FAILED: N=%d\n", n_final))
    overall_results[[disease]] <- data.table(
      disease = disease, design = "landmark",
      HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = n_final, deaths = n_deaths, direction = "ERROR"
    )
    rm(pmap, landmark, cox_pop); gc()
    next
  }

  s <- summary(fit)
  hr <- s$conf.int["E_z", ]
  pval <- s$coefficients["E_z", "Pr(>|z|)"]
  direction <- ifelse(hr[1] < 1, "PASS", "FAIL")

  overall_results[[disease]] <- data.table(
    disease = disease, design = "landmark",
    HR = hr[1], HR_lo = hr[3], HR_hi = hr[4],
    pval = pval, N = n_final, deaths = n_deaths,
    direction = direction,
    n_cohort = n_patients, n_eligible = n_eligible_pts,
    n_excl_riy = n_excluded_riy, n_excl_died_before = n_died_before,
    n_excl_left_before = n_left_before,
    median_fu = median(cox_pop$followup_time),
    median_obs_to_landmark = median(landmark[kojin_id %in% cox_pop$kojin_id,
                                             parse_ym(landmark_ym) - parse_ym(observable_start_ym)])
  )

  cat(sprintf("  N=%d, deaths=%d, HR=%.4f [%.4f, %.4f] → %s\n",
              n_final, n_deaths, hr[1], hr[3], hr[4], direction))

  # --- Per-basin Cox ---
  basins <- sort(unique(cox_pop$basin_id))
  for (b in basins) {
    sub <- cox_pop[basin_id == b]
    n_b <- nrow(sub)
    d_b <- sum(sub$event)

    if (d_b < 10) {
      basin_results[[paste0(disease, "_", b)]] <- data.table(
        disease = disease, basin_id = b, N = n_b, deaths = d_b,
        HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
        pval = NA_real_, direction = "SKIP"
      )
      next
    }

    sub[, E_z_b := scale(energy)[,1]]
    fit_b <- tryCatch(
      coxph(Surv(followup_time, event) ~ E_z_b + age_at_landmark + factor(sex_code),
            data = sub),
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      basin_results[[paste0(disease, "_", b)]] <- data.table(
        disease = disease, basin_id = b, N = n_b, deaths = d_b,
        HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
        pval = NA_real_, direction = "ERROR"
      )
      next
    }

    s_b <- summary(fit_b)
    hr_b <- s_b$conf.int["E_z_b", ]
    pval_b <- s_b$coefficients["E_z_b", "Pr(>|z|)"]
    dir_b <- ifelse(hr_b[1] < 1, "PASS", "FAIL")

    basin_results[[paste0(disease, "_", b)]] <- data.table(
      disease = disease, basin_id = b, N = n_b, deaths = d_b,
      HR = hr_b[1], HR_lo = hr_b[3], HR_hi = hr_b[4],
      pval = pval_b, direction = dir_b
    )
  }

  rm(pmap, landmark, cox_pop)
  gc()

  t1 <- proc.time()["elapsed"]
  cat(sprintf("  done in %.0fs\n", t1 - t0))
}

# ============================================================
# STEP 5: Save results
# ============================================================
cat("\n=== Saving results ===\n")

overall_dt <- rbindlist(overall_results, fill = TRUE)
basin_dt <- rbindlist(basin_results, fill = TRUE)

fwrite(overall_dt, file.path(SCRATCH, "landmark_all_diseases_overall.tsv"), sep = "\t")
fwrite(basin_dt, file.path(SCRATCH, "landmark_all_diseases_per_basin.tsv"), sep = "\t")

cat(sprintf("\nOverall summary:\n"))
cat(sprintf("  PASS: %d diseases\n", sum(overall_dt$direction == "PASS", na.rm = TRUE)))
cat(sprintf("  FAIL: %d diseases\n", sum(overall_dt$direction == "FAIL", na.rm = TRUE)))
cat(sprintf("  SKIP: %d diseases\n", sum(overall_dt$direction == "SKIP", na.rm = TRUE)))
cat(sprintf("  ERROR: %d diseases\n", sum(overall_dt$direction == "ERROR", na.rm = TRUE)))

cat(sprintf("\nPer-basin summary:\n"))
cat(sprintf("  PASS basins: %d\n", sum(basin_dt$direction == "PASS", na.rm = TRUE)))
cat(sprintf("  FAIL basins: %d\n", sum(basin_dt$direction == "FAIL", na.rm = TRUE)))

# ============================================================
# STEP 6: Write pipeline-compatible output files
# ============================================================
cat("\n=== Writing pipeline-compatible output ===\n")

for (di in seq_len(nrow(overall_dt))) {
  row <- overall_dt[di]
  disease <- row$disease
  if (row$direction %in% c("SKIP", "ERROR")) next

  out_dir <- file.path("output", disease, "all_comorbid", "cox", "7")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Per-basin data for this disease
  bdt <- basin_dt[disease == row$disease & !is.na(HR)]
  n_basins <- nrow(bdt)
  largest_frac <- if (n_basins > 0) max(bdt$N) / row$N else NA_real_

  # Overall sanity check (pipeline format)
  sc <- data.table(
    disease_family = disease,
    dim = 7L,
    outcome = "mortality",
    HR = row$HR,
    HR_lower = row$HR_lo,
    HR_upper = row$HR_hi,
    pvalue = row$pval,
    direction_ok = row$direction == "PASS",
    n_patients = row$N,
    n_basins = n_basins,
    largest_basin_frac = largest_frac
  )
  fwrite(sc, file.path(out_dir, "energy_sanity_check.tsv"), sep = "\t")

  # Per-basin sanity check (pipeline format)
  if (n_basins > 0) {
    pb <- data.table(
      disease_family = disease,
      basin_id = bdt$basin_id,
      outcome = "mortality",
      HR = bdt$HR,
      HR_lower = bdt$HR_lo,
      HR_upper = bdt$HR_hi,
      pvalue = bdt$pval,
      direction_ok = bdt$direction == "PASS",
      n_patients = bdt$N,
      n_deaths = bdt$deaths,
      death_rate = bdt$deaths / bdt$N,
      note = ""
    )
    fwrite(pb, file.path(out_dir, "energy_sanity_check_per_basin.tsv"), sep = "\t")
  }
}

# Atlas summary (pipeline format)
atlas <- overall_dt[direction %in% c("PASS", "FAIL")]
atlas_out <- data.table(
  disease_family = atlas$disease,
  dim = 7L,
  outcome = "mortality",
  HR = atlas$HR,
  HR_lower = atlas$HR_lo,
  HR_upper = atlas$HR_hi,
  pvalue = atlas$pval,
  direction_ok = atlas$direction == "PASS",
  n_patients = atlas$N,
  n_basins = sapply(atlas$disease, function(d) {
    sum(!is.na(basin_dt[disease == d, HR]))
  }),
  largest_basin_frac = sapply(atlas$disease, function(d) {
    bdt <- basin_dt[disease == d & !is.na(HR)]
    a_n <- overall_dt[disease == d, N]
    if (nrow(bdt) > 0 && length(a_n) > 0) max(bdt$N) / a_n else NA_real_
  })
)
fwrite(atlas_out, "summary/atlas_summary.tsv", sep = "\t")
cat(sprintf("  Wrote atlas_summary.tsv: %d diseases\n", nrow(atlas_out)))

cat("\n=== ALL DISEASES COMPLETE ===\n")
