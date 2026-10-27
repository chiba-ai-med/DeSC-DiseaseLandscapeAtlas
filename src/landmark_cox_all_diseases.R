## Unified landmark Cox analysis for all 51 diseases
## Builds shared cohort table, runs mortality + hospitalization Cox, generates KM
##
## Usage: Rscript src/landmark_cox_all_diseases.R [window_lookup.bin]

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

bin_arg <- commandArgs(trailingOnly = TRUE)
BIN_PATH <- if (length(bin_arg) > 0) bin_arg[1] else
  "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad/window_lookup.bin"

DATA_START <- "2014/04"
DATA_END   <- "2022/02"
CASE_DISEASES <- c("dementia", "thyroid_disorder", "anemia", "tuberculosis")

# ============================================================
# Utility functions
# ============================================================
parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}

decode_ym <- function(enc) {
  yr <- enc %/% 12L
  mo <- enc %% 12L + 1L
  sprintf("%04d/%02d", yr, mo)
}

# Round-trip test (encode = year*12 + (month-1), decode = %/%12 and %%12+1)
encode_ym <- function(s) {
  as.integer(substr(s, 1, 4)) * 12L + as.integer(substr(s, 6, 7)) - 1L
}

cat("=== Round-trip encode/decode test ===\n")
for (y in c(2014L, 2020L, 2022L)) {
  for (m in 1:12) {
    s <- sprintf("%04d/%02d", y, m)
    enc <- encode_ym(s)
    dec <- decode_ym(enc)
    stopifnot(dec == s)
  }
}
cat("  PASSED: all 36 year-month combinations\n")

# ============================================================
# STEP 1: Load compact binary window lookup
# ============================================================
cat("\n=== Loading compact window lookup ===\n")
stopifnot(file.exists(BIN_PATH))

N_ROWS <- 347764079L
con <- file(BIN_PATH, "rb")
raw <- readBin(con, "integer", n = N_ROWS * 2L, size = 2L,
               signed = FALSE, endian = "little")
close(con)
stopifnot(length(raw) == N_ROWS * 2L)

win_start_enc <- raw[seq(1L, length(raw), 2L)]
win_end_enc   <- raw[seq(2L, length(raw), 2L)]
rm(raw); gc()
cat(sprintf("  Loaded %d window pairs\n", length(win_start_enc)))

# Verify row_id alignment: spot-check a few known rows
# Row 1 should decode to 2014/04..2014/09 (from the TSV)
r1_start <- decode_ym(win_start_enc[1])
r1_end   <- decode_ym(win_end_enc[1])
cat(sprintf("  Row 1: %s..%s (expected 2014/04..2014/09)\n", r1_start, r1_end))
stopifnot(r1_start == "2014/04" && r1_end == "2014/09")

# ============================================================
# STEP 2: Load tekiyo
# ============================================================
cat("\n=== Loading tekiyo ===\n")
tek <- fread("data/tekiyo.csv",
             select = c("kojin_id", "birth_ym", "sex_code",
                        "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
tek[, shibou_flg_riyouka := as.integer(shibou_flg_riyouka)]
# For riyouka=1: shibou_flg "1" → 1 (dead), anything else → 0 (alive/censored)
# For riyouka=0/NA: leave shibou_flg as-is (not used in mortality Cox)
tek[shibou_flg_riyouka == 1L,
    event_death := fifelse(!is.na(shibou_flg) & shibou_flg == "1", 1L, 0L)]
cat(sprintf("  %d patients, %d with riyouka=1, %d deaths among riyouka=1\n",
            nrow(tek),
            sum(tek$shibou_flg_riyouka == 1L, na.rm = TRUE),
            sum(tek$event_death == 1L, na.rm = TRUE)))

# ============================================================
# STEP 3: Disease list
# ============================================================
diseases <- fread("config/disease_families.tsv", sep = "\t")
disease_list <- diseases[family_id != "hf_cancer", family_id]
cat(sprintf("\n=== Processing %d diseases ===\n\n", length(disease_list)))

# ============================================================
# STEP 4: Process each disease
# ============================================================
overall_mort <- list()
overall_hosp <- list()
basin_mort <- list()
basin_hosp <- list()

for (di in seq_along(disease_list)) {
  disease <- disease_list[di]
  base_dir <- file.path("output", disease, "all_comorbid")
  cat(sprintf("[%d/%d] %s\n", di, length(disease_list), disease))
  t0 <- proc.time()["elapsed"]

  # --- Load person_map ---
  pmap_path <- file.path(base_dir, "cohort", "person_map.tsv")
  pmap <- fread(pmap_path)
  pmap[, kojin_id := as.character(kojin_id)]

  # --- Window lookup via binary index ---
  valid_idx <- pmap$old_row_id >= 1L & pmap$old_row_id <= N_ROWS
  pmap[valid_idx, `:=`(
    ws_enc = win_start_enc[old_row_id],
    we_enc = win_end_enc[old_row_id]
  )]
  pmap <- pmap[!is.na(ws_enc) & ws_enc > 0L & !is.na(we_enc) & we_enc > 0L]
  pmap[, window_start := decode_ym(ws_enc)]
  pmap[, window_end   := decode_ym(we_enc)]
  pmap[, c("ws_enc", "we_enc") := NULL]

  # --- Landscape: bit-key → basin_id + energy ---
  ls_dir <- file.path(base_dir, "landscaper", "7")
  sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))

  scores_mm_path <- file.path(base_dir, "pca", "7", "Scores.mm")
  if (!file.exists(scores_mm_path)) {
    cat("  SKIP: no Scores.mm\n"); next
  }

  hdr <- readLines(scores_mm_path, n = 2)
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

  key_dt[, pattern_id := as.integer(key2pid[as.character(key)])]
  setnames(key_dt, "i", "new_row_id")
  key_dt <- merge(key_dt, sbmap[, .(pattern_id = state_id, basin_id, energy)],
                  by = "pattern_id", all.x = TRUE)

  pmap <- merge(pmap, key_dt[, .(new_row_id, basin_id, energy)],
                by = "new_row_id", all.x = TRUE)
  rm(key_dt)

  # All-negative state (key=0, absent from Scores.mm)
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

  # --- Eligibility: basin_id > 0, window within insurance + data period ---
  pmap[, eligible := (
    !is.na(basin_id) & basin_id > 0L &
    window_start >= observable_start_ym &
    window_end <= observable_end_ym &
    window_start >= DATA_START &
    window_end <= DATA_END
  )]

  eligible_obs <- pmap[eligible == TRUE]
  if (nrow(eligible_obs) == 0) {
    cat("  SKIP: no eligible observations\n")
    overall_mort[[disease]] <- data.table(
      disease = disease, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = 0L, deaths = 0L, censored = 0L, direction = "SKIP")
    rm(pmap, eligible_obs); gc(); next
  }

  # --- First eligible window by actual calendar month ---
  setorder(eligible_obs, kojin_id, window_start, old_row_id)
  landmark <- eligible_obs[!duplicated(kojin_id)]
  rm(eligible_obs)

  # --- Follow-up: from window_end to observable_end_ym ---
  landmark[, fu_time := parse_ym(observable_end_ym) - parse_ym(window_end)]
  landmark <- landmark[fu_time > 0]

  # Age at landmark
  landmark[, age_at_lm := parse_ym(window_start) - parse_ym(birth_ym)]

  # --- Load disease-specific hospitalization ---
  hosp_path <- file.path(base_dir, "cohort", "person_hospitalization_specific.tsv")
  if (file.exists(hosp_path)) {
    hosp <- fread(hosp_path)
    hosp[, kojin_id := as.character(kojin_id)]
    hosp[, admission_ym := substr(first_admission_ymd, 1, 7)]
    landmark <- merge(landmark, hosp[, .(kojin_id, admission_ym)],
                      by = "kojin_id", all.x = TRUE)
    rm(hosp)
  } else {
    landmark[, admission_ym := NA_character_]
  }

  # Hospitalization event: admission after window_end
  landmark[, event_hosp := fifelse(
    !is.na(admission_ym) & admission_ym > window_end, 1L, 0L)]
  # Hospitalization follow-up: to admission if post-landmark, else to obs_end
  landmark[, fu_hosp := fifelse(
    event_hosp == 1L,
    parse_ym(admission_ym) - parse_ym(window_end),
    fu_time)]
  landmark[fu_hosp <= 0, `:=`(event_hosp = 0L, fu_hosp = fu_time)]

  n_all <- nrow(landmark)

  # --- Save landmark cohort table ---
  cohort_dir <- file.path(base_dir, "landmark")
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
  cohort_out <- landmark[, .(
    kojin_id, window_start, window_end, basin_id, energy,
    birth_ym, sex_code, observable_start_ym, observable_end_ym,
    shibou_flg_riyouka, event_death, fu_time,
    admission_ym, event_hosp, fu_hosp, age_at_lm)]
  fwrite(cohort_out, file.path(cohort_dir, "landmark_cohort.tsv"), sep = "\t")
  rm(cohort_out)

  # ==========================================================
  # MORTALITY COX (riyouka=1 subset)
  # ==========================================================
  mort_pop <- landmark[shibou_flg_riyouka == 1L]
  n_mort <- nrow(mort_pop)
  n_deaths <- sum(mort_pop$event_death, na.rm = TRUE)
  n_cens_mort <- n_mort - n_deaths

  if (n_mort < 100 || n_deaths < 20) {
    cat(sprintf("  Mortality SKIP: N=%d, deaths=%d\n", n_mort, n_deaths))
    overall_mort[[disease]] <- data.table(
      disease = disease, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = n_mort, deaths = n_deaths, censored = n_cens_mort,
      direction = "SKIP")
  } else {
    mort_pop[, E_z := scale(energy)[,1]]

    fit <- tryCatch(
      coxph(Surv(fu_time, event_death) ~ E_z + age_at_lm + factor(sex_code),
            data = mort_pop),
      error = function(e) NULL)

    if (is.null(fit)) {
      overall_mort[[disease]] <- data.table(
        disease = disease, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
        pval = NA_real_, N = n_mort, deaths = n_deaths, censored = n_cens_mort,
        direction = "ERROR")
    } else {
      s <- summary(fit)
      hr <- s$conf.int["E_z", ]
      pval <- s$coefficients["E_z", "Pr(>|z|)"]
      dir_ <- ifelse(hr[1] < 1, "PASS", "FAIL")

      overall_mort[[disease]] <- data.table(
        disease = disease, HR = hr[1], HR_lo = hr[3], HR_hi = hr[4],
        pval = pval, N = n_mort, deaths = n_deaths, censored = n_cens_mort,
        direction = dir_,
        median_fu = median(mort_pop$fu_time))

      cat(sprintf("  Mortality: N=%d, deaths=%d, cens=%d, HR=%.4f [%.4f,%.4f] → %s\n",
                  n_mort, n_deaths, n_cens_mort, hr[1], hr[3], hr[4], dir_))

      # Per-basin mortality Cox
      for (b in sort(unique(mort_pop$basin_id))) {
        sub <- mort_pop[basin_id == b]
        d_b <- sum(sub$event_death)
        if (d_b < 10) {
          basin_mort[[paste0(disease,"_",b)]] <- data.table(
            disease=disease, basin_id=b, N=nrow(sub), deaths=d_b,
            censored=nrow(sub)-d_b,
            HR=NA_real_, HR_lo=NA_real_, HR_hi=NA_real_,
            pval=NA_real_, direction="SKIP")
          next
        }
        sub[, E_z_b := scale(energy)[,1]]
        fit_b <- tryCatch(
          coxph(Surv(fu_time, event_death) ~ E_z_b + age_at_lm + factor(sex_code),
                data = sub),
          error = function(e) NULL)
        if (is.null(fit_b)) {
          basin_mort[[paste0(disease,"_",b)]] <- data.table(
            disease=disease, basin_id=b, N=nrow(sub), deaths=d_b,
            censored=nrow(sub)-d_b,
            HR=NA_real_, HR_lo=NA_real_, HR_hi=NA_real_,
            pval=NA_real_, direction="ERROR")
          next
        }
        s_b <- summary(fit_b)
        hr_b <- s_b$conf.int["E_z_b",]
        pv_b <- s_b$coefficients["E_z_b","Pr(>|z|)"]
        basin_mort[[paste0(disease,"_",b)]] <- data.table(
          disease=disease, basin_id=b, N=nrow(sub), deaths=d_b,
          censored=nrow(sub)-d_b,
          HR=hr_b[1], HR_lo=hr_b[3], HR_hi=hr_b[4],
          pval=pv_b, direction=ifelse(hr_b[1]<1,"PASS","FAIL"))
      }
    }
  }
  rm(mort_pop)

  # ==========================================================
  # HOSPITALIZATION COX (all eligible patients)
  # ==========================================================
  n_hosp_events <- sum(landmark$event_hosp)

  if (n_all < 100 || n_hosp_events < 20) {
    cat(sprintf("  Hosp SKIP: N=%d, events=%d\n", n_all, n_hosp_events))
    overall_hosp[[disease]] <- data.table(
      disease = disease, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      pval = NA_real_, N = n_all, events = n_hosp_events,
      censored = n_all - n_hosp_events, direction = "SKIP")
  } else {
    landmark[, E_z_h := scale(energy)[,1]]
    fit_h <- tryCatch(
      coxph(Surv(fu_hosp, event_hosp) ~ E_z_h + age_at_lm + factor(sex_code),
            data = landmark),
      error = function(e) NULL)

    if (is.null(fit_h)) {
      overall_hosp[[disease]] <- data.table(
        disease = disease, HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
        pval = NA_real_, N = n_all, events = n_hosp_events,
        censored = n_all - n_hosp_events, direction = "ERROR")
    } else {
      s_h <- summary(fit_h)
      hr_h <- s_h$conf.int["E_z_h",]
      pv_h <- s_h$coefficients["E_z_h","Pr(>|z|)"]
      dir_h <- ifelse(hr_h[1] < 1, "PASS", "FAIL")

      overall_hosp[[disease]] <- data.table(
        disease = disease, HR = hr_h[1], HR_lo = hr_h[3], HR_hi = hr_h[4],
        pval = pv_h, N = n_all, events = n_hosp_events,
        censored = n_all - n_hosp_events, direction = dir_h)

      cat(sprintf("  Hosp: N=%d, events=%d, HR=%.4f → %s\n",
                  n_all, n_hosp_events, hr_h[1], dir_h))

      # Per-basin hospitalization Cox
      for (b in sort(unique(landmark$basin_id))) {
        sub <- landmark[basin_id == b]
        ev_b <- sum(sub$event_hosp)
        if (ev_b < 10) {
          basin_hosp[[paste0(disease,"_",b)]] <- data.table(
            disease=disease, basin_id=b, N=nrow(sub), events=ev_b,
            censored=nrow(sub)-ev_b,
            HR=NA_real_, HR_lo=NA_real_, HR_hi=NA_real_,
            pval=NA_real_, direction="SKIP")
          next
        }
        sub[, E_z_hb := scale(energy)[,1]]
        fit_hb <- tryCatch(
          coxph(Surv(fu_hosp, event_hosp) ~ E_z_hb + age_at_lm + factor(sex_code),
                data = sub),
          error = function(e) NULL)
        if (is.null(fit_hb)) {
          basin_hosp[[paste0(disease,"_",b)]] <- data.table(
            disease=disease, basin_id=b, N=nrow(sub), events=ev_b,
            censored=nrow(sub)-ev_b,
            HR=NA_real_, HR_lo=NA_real_, HR_hi=NA_real_,
            pval=NA_real_, direction="ERROR")
          next
        }
        s_hb <- summary(fit_hb)
        hr_hb <- s_hb$conf.int["E_z_hb",]
        pv_hb <- s_hb$coefficients["E_z_hb","Pr(>|z|)"]
        basin_hosp[[paste0(disease,"_",b)]] <- data.table(
          disease=disease, basin_id=b, N=nrow(sub), events=ev_b,
          censored=nrow(sub)-ev_b,
          HR=hr_hb[1], HR_lo=hr_hb[3], HR_hi=hr_hb[4],
          pval=pv_hb, direction=ifelse(hr_hb[1]<1,"PASS","FAIL"))
      }
    }
  }

  # ==========================================================
  # KM CURVES (case study diseases only)
  # ==========================================================
  if (disease %in% CASE_DISEASES) {
    cat("  Generating KM curves...\n")

    # --- Mortality KM (riyouka=1 subset) ---
    km_pop <- landmark[shibou_flg_riyouka == 1L]
    km_pop[, basin_label := paste0("Basin ", basin_id)]

    survfit_to_dt <- function(sf) {
      strata_names <- if (!is.null(sf$strata)) {
        rep(names(sf$strata), sf$strata)
      } else {
        rep("All", length(sf$time))
      }
      data.table(time = sf$time, surv = sf$surv,
                 upper = sf$upper, lower = sf$lower,
                 strata = strata_names)
    }

    plot_dir <- file.path("plot", disease)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    # Inter-basin KM
    if (length(unique(km_pop$basin_id)) > 1) {
      fit_km <- survfit(Surv(fu_time, event_death) ~ basin_label, data = km_pop)
      km_dt <- survfit_to_dt(fit_km)
      km_dt[, strata := gsub("basin_label=", "", strata)]
      pv <- tryCatch({
        sd <- survdiff(Surv(fu_time, event_death) ~ basin_label, data = km_pop)
        1 - pchisq(sd$chisq, length(sd$n) - 1)
      }, error = function(e) NA)
      pv_label <- if (is.na(pv)) "" else if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)

      g1 <- ggplot(km_dt, aes(x = time, y = surv, color = strata)) +
        geom_step(linewidth = 0.8) +
        labs(title = sprintf("%s — Landmark KM by basin (mortality)",
                             gsub("_", " ", disease)),
             subtitle = pv_label,
             x = "Years from landmark", y = "Survival probability",
             color = "Basin") +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom") +
        coord_cartesian(xlim = c(0, min(8, max(km_dt$time))))
      ggsave(file.path(plot_dir, "km_inter_basin_landmark.png"),
             plot = g1, width = 8, height = 6, dpi = 200)
      cat(sprintf("  Written: %s/km_inter_basin_landmark.png\n", plot_dir))
    }

    # Per-basin energy tertile KM
    plots <- list()
    for (b in sort(unique(km_pop$basin_id))) {
      sub <- km_pop[basin_id == b]
      if (nrow(sub) < 50) next
      brks <- quantile(sub$energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
      if (length(unique(brks)) < 4) next
      sub[, e_tert := cut(energy, breaks = brks,
        labels = c("Low E", "Mid E", "High E"), include.lowest = TRUE)]
      fit_b <- survfit(Surv(fu_time, event_death) ~ e_tert, data = sub)
      km_b <- survfit_to_dt(fit_b)
      km_b[, strata := gsub("e_tert=", "", strata)]
      km_b[, basin := sprintf("Basin %d (n=%s)", b, format(nrow(sub), big.mark=","))]
      plots[[length(plots) + 1]] <- km_b
    }
    if (length(plots) > 0) {
      all_km <- rbindlist(plots)
      g2 <- ggplot(all_km, aes(x = time, y = surv, color = strata)) +
        geom_step(linewidth = 0.7) +
        facet_wrap(~ basin, nrow = 1) +
        scale_color_manual(values = c("Low E" = "#2166AC", "Mid E" = "#BDBDBD",
                                       "High E" = "#B2182B")) +
        labs(x = "Years from landmark", y = "Survival", color = "Energy tertile",
             title = sprintf("%s — Within-basin energy tertile KM",
                             gsub("_", " ", disease))) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom") +
        coord_cartesian(xlim = c(0, min(8, max(all_km$time))))
      ggsave(file.path(plot_dir, "km_per_basin_landmark.png"), plot = g2,
             width = 4 * min(length(plots), 4), height = 5, dpi = 200)
      cat(sprintf("  Written: %s/km_per_basin_landmark.png\n", plot_dir))
    }
    rm(km_pop)
  }

  rm(pmap, landmark); gc()
  t1 <- proc.time()["elapsed"]
  cat(sprintf("  done in %.0fs\n", t1 - t0))
}

# ============================================================
# STEP 5: Save results
# ============================================================
cat("\n=== Saving results ===\n")
mort_dt <- rbindlist(overall_mort, fill = TRUE)
hosp_dt <- rbindlist(overall_hosp, fill = TRUE)
basin_mort_dt <- rbindlist(basin_mort, fill = TRUE)
basin_hosp_dt <- rbindlist(basin_hosp, fill = TRUE)

fwrite(mort_dt, "summary/landmark_all_diseases_overall.tsv", sep = "\t")
fwrite(hosp_dt, "summary/landmark_all_diseases_hosp_overall.tsv", sep = "\t")
fwrite(basin_mort_dt, "summary/landmark_all_diseases_per_basin.tsv", sep = "\t")
fwrite(basin_hosp_dt, "summary/landmark_all_diseases_hosp_per_basin.tsv", sep = "\t")

cat(sprintf("\nMortality: PASS=%d, FAIL=%d, SKIP=%d\n",
    sum(mort_dt$direction == "PASS", na.rm = TRUE),
    sum(mort_dt$direction == "FAIL", na.rm = TRUE),
    sum(mort_dt$direction == "SKIP", na.rm = TRUE)))
cat(sprintf("Hospitalization: PASS=%d, FAIL=%d, SKIP=%d\n",
    sum(hosp_dt$direction == "PASS", na.rm = TRUE),
    sum(hosp_dt$direction == "FAIL", na.rm = TRUE),
    sum(hosp_dt$direction == "SKIP", na.rm = TRUE)))

# ============================================================
# STEP 6: Post-process — atlas_summary.tsv
# ============================================================
cat("\n=== Building atlas_summary.tsv ===\n")

# Per-basin majority criterion for mortality (basin_id > 0 only, already guaranteed)
pb_mort <- basin_mort_dt[direction %in% c("PASS", "FAIL"),
  .(n_pass = sum(direction == "PASS"), n_basins = .N), by = disease]
pb_mort[, majority_pass := n_pass >= ceiling(n_basins / 2)]

# Per-basin majority criterion for hospitalization
pb_hosp <- basin_hosp_dt[direction %in% c("PASS", "FAIL"),
  .(n_pass_h = sum(direction == "PASS"), n_basins_h = .N), by = disease]
pb_hosp[, majority_pass_h := n_pass_h >= ceiling(n_basins_h / 2)]

# Build atlas
atlas <- mort_dt[direction %in% c("PASS", "FAIL")]
atlas <- merge(atlas, pb_mort[, .(disease, n_basins, majority_pass)],
               by = "disease", all.x = TRUE)
atlas <- merge(atlas, pb_hosp[, .(disease, majority_pass_h)],
               by = "disease", all.x = TRUE)

# Largest basin fraction
atlas[, largest_basin_frac := sapply(disease, function(d) {
  bdt <- basin_mort_dt[disease == d & !is.na(HR)]
  a_n <- mort_dt[disease == d, N]
  if (nrow(bdt) > 0 && length(a_n) > 0 && a_n > 0) max(bdt$N) / a_n else NA_real_
})]

atlas_out <- data.table(
  disease_family = atlas$disease,
  dim = 7L,
  outcome = "mortality",
  HR = atlas$HR,
  HR_lower = atlas$HR_lo,
  HR_upper = atlas$HR_hi,
  pvalue = atlas$pval,
  n_patients = atlas$N,
  n_basins = atlas$n_basins,
  largest_basin_frac = atlas$largest_basin_frac,
  direction_ok = atlas$majority_pass
)

fwrite(atlas_out, "summary/atlas_summary.tsv", sep = "\t")
cat(sprintf("  atlas_summary.tsv: %d diseases\n", nrow(atlas_out)))
cat(sprintf("  Mortality per-basin majority: PASS=%d, FAIL=%d\n",
    sum(atlas_out$direction_ok == TRUE, na.rm = TRUE),
    sum(atlas_out$direction_ok == FALSE, na.rm = TRUE)))

# Hospitalization per-basin majority summary
n_hosp_pass <- sum(pb_hosp$majority_pass_h, na.rm = TRUE)
n_hosp_fail <- sum(!pb_hosp$majority_pass_h, na.rm = TRUE)
cat(sprintf("  Hospitalization per-basin majority: PASS=%d, FAIL=%d\n",
    n_hosp_pass, n_hosp_fail))

# ============================================================
# STEP 7: Atlas heatmap
# ============================================================
cat("\n=== Generating atlas heatmap ===\n")
sm <- fread("summary/atlas_summary.tsv")
sm[, direction_label := fifelse(direction_ok, "PASS", "FAIL")]
sm[, hr_label := sprintf("%.3f", HR)]
sm[, disease_family := factor(disease_family,
     levels = rev(sort(unique(disease_family))))]

g_atlas <- ggplot(sm, aes(x = outcome, y = disease_family, fill = direction_label)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = hr_label), size = 3) +
  scale_fill_manual(
    values = c(PASS = "#4CAF50", FAIL = "#F44336"),
    name = "Per-basin majority"
  ) +
  labs(title = "Energy Cox Sanity Check Atlas (Landmark Design)",
       subtitle = "Overall HR; color = per-basin majority (HR<1 in ≥half basins, basin_id>0)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

ggsave("plot/atlas_heatmap.png", g_atlas, width = 8, height = 16, dpi = 200)
cat("  Written: plot/atlas_heatmap.png\n")

cat("\n=== ALL COMPLETE ===\n")
