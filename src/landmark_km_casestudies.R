## Landmark KM curves for case study diseases
## Mirrors the data loading approach of landmark_cox_all_diseases.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

BASE <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
SCRATCH <- "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad"
setwd(BASE)

CASE_DISEASES <- c("dementia", "thyroid_disorder", "anemia", "tuberculosis")

parse_ym <- function(ym_str) {
  parts <- strsplit(ym_str, "/")
  sapply(parts, function(p) as.numeric(p[1]) + (as.numeric(p[2]) - 1) / 12)
}
decode_ym <- function(enc) {
  yr <- enc %/% 12L
  mo <- enc %% 12L
  sprintf("%04d/%02d", yr, mo)
}

## ─── Load shared data once ───
cat("Loading compact binary window lookup...\n")
con <- file(file.path(SCRATCH, "window_lookup.bin"), "rb")
raw <- readBin(con, "integer", n = 347764079L * 2L, size = 2L,
               signed = FALSE, endian = "little")
close(con)
win_start_enc <- raw[seq(1L, length(raw), 2L)]
win_end_enc   <- raw[seq(2L, length(raw), 2L)]
rm(raw); gc()
cat(sprintf("  Loaded %d window pairs\n", length(win_start_enc)))

cat("Loading tekiyo...\n")
tek <- fread("data/tekiyo.csv",
             select = c("kojin_id", "birth_ym", "sex_code",
                        "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
tek[, shibou_flg := fifelse(shibou_flg == "1", 1L, 0L)]
tek[, shibou_flg_riyouka := as.integer(shibou_flg_riyouka)]
cat(sprintf("  %d patients\n", nrow(tek)))

## ─── Process each case study ───
for (disease in CASE_DISEASES) {
  cat(sprintf("\n========== %s ==========\n", disease))
  t0 <- Sys.time()

  ## Load person_map
  pm_path <- file.path("output", disease, "all_comorbid", "cohort", "person_map.tsv")
  pmap <- fread(pm_path, select = c("kojin_id", "old_row_id"))
  pmap[, kojin_id := as.character(kojin_id)]
  cat(sprintf("  person_map: %d rows\n", nrow(pmap)))

  ## Decode windows
  pmap[, ws_enc := win_start_enc[old_row_id]]
  pmap[, we_enc := win_end_enc[old_row_id]]
  pmap <- pmap[ws_enc > 0L & we_enc > 0L]
  pmap[, window_start := decode_ym(ws_enc)]
  pmap[, window_end   := decode_ym(we_enc)]
  pmap[, c("ws_enc", "we_enc") := NULL]

  ## Load Scores.mm
  scores_mm_path <- file.path("output", disease, "all_comorbid", "pca", "7", "Scores.mm")
  hdr <- readLines(scores_mm_path, n = 2)
  dims <- as.integer(strsplit(trimws(hdr[2]), " +")[[1]])
  ij_dt <- fread(scores_mm_path, skip = 2, header = FALSE,
                 col.names = c("i", "j"), select = 1:2)
  pos_cols <- ij_dt[, .(pos_cols = list(j)), by = i]
  setnames(pos_cols, "i", "obs_row")
  rm(ij_dt); gc()

  ## Row mapping
  row_ids_path <- file.path("output", disease, "all_comorbid", "cohort", "row_ids.txt")
  row_map <- fread(row_ids_path, header = FALSE, col.names = c("old_row_id", "new_row_id"))
  pmap <- merge(pmap, row_map, by = "old_row_id", all.x = TRUE)
  pmap <- pmap[!is.na(new_row_id)]
  rm(row_map)

  ## Merge PCA positive columns & bit encode
  pmap <- merge(pmap, pos_cols, by.x = "new_row_id", by.y = "obs_row", all.x = TRUE)
  rm(pos_cols)
  pmap[, bitkey := sapply(pos_cols, function(cols) {
    if (is.null(cols) || length(cols) == 0) return(0L)
    sum(2L^(cols - 1L))
  })]
  pmap[, pos_cols := NULL]

  ## Load landscape via state_basin_map.tsv
  ls_dir <- file.path("output", disease, "all_comorbid", "landscaper", "7")
  sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))
  allstates <- fread(file.path(ls_dir, "Allstates.tsv"))

  as_mat <- as.matrix(allstates)
  keys <- apply(as_mat, 1, function(row) {
    pos <- which(row == 1)
    if (length(pos) == 0) 0L else sum(2L^(pos - 1L))
  })
  key_dt <- data.table(bitkey = keys, state_id = seq_along(keys))
  key_dt <- merge(key_dt, sbmap[, .(pattern_id = state_id, basin_id, energy)],
                  by.x = "state_id", by.y = "pattern_id", all.x = TRUE)
  key_dt[is.na(basin_id), basin_id := 0L]

  ## All-negative state
  all_neg_pid <- key_dt[bitkey == 0, state_id]
  if (length(all_neg_pid) > 0) {
    neg_info <- sbmap[state_id == all_neg_pid]
    if (nrow(neg_info) > 0) {
      key_dt[bitkey == 0, `:=`(basin_id = neg_info$basin_id[1],
                                energy = neg_info$energy[1])]
    }
  }

  pmap <- merge(pmap, key_dt[, .(bitkey, energy, basin_id)], by = "bitkey", all.x = TRUE)
  pmap[is.na(basin_id), basin_id := 0L]
  rm(key_dt)

  ## Merge tekiyo
  pmap <- merge(pmap, tek, by = "kojin_id", all.x = TRUE)

  ## Eligibility filter
  pmap <- pmap[shibou_flg_riyouka == 1L]
  pmap[, obs_start_yr := parse_ym(observable_start_ym)]
  pmap[, obs_end_yr   := parse_ym(observable_end_ym)]
  pmap[, win_start_yr := parse_ym(window_start)]
  pmap[, win_end_yr   := parse_ym(window_end)]

  pmap <- pmap[win_start_yr >= obs_start_yr &
               win_end_yr <= obs_end_yr &
               window_start >= "2014/04" &
               window_end <= "2022/02" &
               basin_id > 0L]

  ## First eligible window per patient
  setorder(pmap, kojin_id, old_row_id)
  first <- pmap[!duplicated(kojin_id)]
  rm(pmap); gc()

  ## Follow-up time
  first[, fu_time := obs_end_yr - win_end_yr]
  first <- first[fu_time > 0]
  first[, birth_y := parse_ym(birth_ym)]
  first[, age := win_start_yr - birth_y]

  cat(sprintf("  Landmark eligible: %d patients\n", nrow(first)))
  cat(sprintf("  Deaths: %d (%.1f%%)\n",
              sum(first$shibou_flg == 1), 100 * mean(first$shibou_flg == 1)))

  first[, basin_label := paste0("Basin ", basin_id)]

  ## Helper: extract survfit as data.table
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

  ## ─── KM: inter-basin ───
  if (length(unique(first$basin_id)) > 1) {
    fit_basin <- survfit(Surv(fu_time, shibou_flg) ~ basin_label, data = first)
    km_dt <- survfit_to_dt(fit_basin)
    km_dt[, strata := gsub("basin_label=", "", strata)]
    pv <- tryCatch({
      sd <- survdiff(Surv(fu_time, shibou_flg) ~ basin_label, data = first)
      1 - pchisq(sd$chisq, length(sd$n) - 1)
    }, error = function(e) NA)
    pv_label <- if (is.na(pv)) "" else if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)

    p <- ggplot(km_dt, aes(x = time, y = surv, color = strata)) +
      geom_step(linewidth = 0.8) +
      labs(title = sprintf("%s — Landmark KM by basin (mortality)",
                           gsub("_", " ", disease)),
           subtitle = pv_label,
           x = "Years from landmark", y = "Survival probability",
           color = "Basin") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom") +
      coord_cartesian(xlim = c(0, min(8, max(km_dt$time))))
    out_path <- file.path("plot", disease, "km_inter_basin_landmark.png")
    ggsave(out_path, plot = p, width = 8, height = 6, dpi = 200)
    cat(sprintf("  Written: %s\n", out_path))
  } else {
    cat("  Single basin — no inter-basin KM\n")
  }

  ## ─── KM: energy tertile per basin ───
  basins_present <- sort(unique(first$basin_id))
  plots <- list()
  for (b in basins_present) {
    sub <- first[basin_id == b]
    if (nrow(sub) < 50) next
    brks <- quantile(sub$energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    if (length(unique(brks)) < 4) next
    sub[, e_tert := cut(energy, breaks = brks,
      labels = c("Low E", "Mid E", "High E"), include.lowest = TRUE)]
    fit_b <- survfit(Surv(fu_time, shibou_flg) ~ e_tert, data = sub)
    km_b <- survfit_to_dt(fit_b)
    km_b[, strata := gsub("e_tert=", "", strata)]
    km_b[, basin := sprintf("Basin %d (n=%s)", b, format(nrow(sub), big.mark=","))]
    plots[[length(plots) + 1]] <- km_b
  }
  if (length(plots) > 0) {
    all_km <- rbindlist(plots)
    p2 <- ggplot(all_km, aes(x = time, y = surv, color = strata)) +
      geom_step(linewidth = 0.7) +
      facet_wrap(~ basin, nrow = 1) +
      scale_color_manual(values = c("Low E" = "#2166AC", "Mid E" = "#BDBDBD", "High E" = "#B2182B")) +
      labs(x = "Years from landmark", y = "Survival", color = "Energy tertile",
           title = sprintf("%s — Within-basin energy tertile KM",
                           gsub("_", " ", disease))) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      coord_cartesian(xlim = c(0, min(8, max(all_km$time))))
    out_path2 <- file.path("plot", disease, "km_per_basin_landmark.png")
    ggsave(out_path2, plot = p2,
           width = 4 * min(length(plots), 4), height = 5, dpi = 200)
    cat(sprintf("  Written: %s\n", out_path2))
  }

  dt <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("  Done in %.0fs\n", dt))
  rm(first); gc()
}

cat("\n=== All case study KM complete ===\n")
