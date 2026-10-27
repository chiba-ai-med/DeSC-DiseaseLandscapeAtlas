## Landmark KM curves from shared cohort tables
## Reads output/{disease}/all_comorbid/landmark/landmark_cohort.tsv
## (produced by landmark_cox_all_diseases.R)
##
## Usage: Rscript src/landmark_km_casestudies.R [disease1 disease2 ...]

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

args <- commandArgs(trailingOnly = TRUE)
CASE_DISEASES <- if (length(args) > 0) args else
  c("dementia", "thyroid_disorder", "anemia", "tuberculosis")

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

for (disease in CASE_DISEASES) {
  cat(sprintf("\n========== %s ==========\n", disease))
  t0 <- Sys.time()

  cohort_path <- file.path("output", disease, "all_comorbid",
                           "landmark", "landmark_cohort.tsv")
  if (!file.exists(cohort_path)) {
    cat("  SKIP: no landmark_cohort.tsv (run landmark_cox_all_diseases.R first)\n")
    next
  }
  dt <- fread(cohort_path)
  cat(sprintf("  Cohort: %d patients, %d basins\n",
              nrow(dt), length(unique(dt$basin_id))))

  # Mortality KM (riyouka=1 subset)
  km_pop <- dt[shibou_flg_riyouka == 1L]
  n_deaths <- sum(km_pop$event_death, na.rm = TRUE)
  cat(sprintf("  Mortality pop: %d patients, %d deaths (%.1f%%)\n",
              nrow(km_pop), n_deaths, 100 * n_deaths / nrow(km_pop)))

  km_pop[, basin_label := paste0("Basin ", basin_id)]
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
  } else {
    cat("  Single basin — no inter-basin KM\n")
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

  dt_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  Done in %.0fs\n", dt_sec))
}

cat("\n=== All case study KM complete ===\n")
