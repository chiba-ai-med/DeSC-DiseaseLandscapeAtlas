## Basin-stratified survival analysis
##
## For each basin: KM curves by energy tertile + Cox HR
##
## Usage:
##   Rscript src/cox_basin_stratified.R <disease> <patient_scores_tsv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease     <- args[1]
scores_path <- args[2]
out_dir     <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Basin-stratified Cox: %s ===\n", disease))

pt <- fread(scores_path)
pt[, death_bin := as.integer(death == "dead")]
pt[, time_years := as.numeric(obs_years)]
pt[, sex_f := factor(sex)]

pt <- pt[is.finite(time_years) & time_years > 0 &
         is.finite(energy) & is.finite(age) & !is.na(basin_id)]

basins <- sort(unique(pt$basin_id))
n_basins <- length(basins)
cat(sprintf("  %d patients, %d basins, %d deaths (%.1f%%)\n",
            nrow(pt), n_basins, sum(pt$death_bin),
            100 * mean(pt$death_bin)))

if (n_basins < 1) {
  cat("No basins, skipping\n")
  quit(save = "no")
}

# --- Basin-stratified Cox results ---
cox_results <- list()

for (b in basins) {
  sub <- pt[basin_id == b]
  n_deaths <- sum(sub$death_bin)
  cat(sprintf("  Basin %d: n=%d, deaths=%d\n", b, nrow(sub), n_deaths))

  if (n_deaths >= 10 && nrow(sub) >= 30) {
    sub[, E_z := scale(energy)[, 1]]
    fit <- coxph(Surv(time_years, death_bin) ~ E_z + age + sex_f, data = sub)
    s <- summary(fit)
    hr <- s$conf.int["E_z", ]
    pv <- coef(s)["E_z", 5]
    cox_results[[length(cox_results) + 1]] <- data.table(
      disease_family = disease, basin_id = b,
      n_patients = nrow(sub), n_deaths = n_deaths,
      HR = hr[1], HR_lower = hr[3], HR_upper = hr[4],
      pvalue = pv, direction_ok = (hr[1] < 1.0)
    )
    cat(sprintf("    HR(E_z) = %.4f [%.4f, %.4f] p=%.2e -> %s\n",
                hr[1], hr[3], hr[4], pv,
                ifelse(hr[1] < 1, "PASS", "FAIL")))
  } else {
    cox_results[[length(cox_results) + 1]] <- data.table(
      disease_family = disease, basin_id = b,
      n_patients = nrow(sub), n_deaths = n_deaths,
      HR = NA_real_, HR_lower = NA_real_, HR_upper = NA_real_,
      pvalue = NA_real_, direction_ok = NA
    )
    cat("    Too few events\n")
  }
}

cox_dt <- rbindlist(cox_results)
fwrite(cox_dt, file.path(out_dir, "cox_basin_stratified.tsv"), sep = "\t")

# --- KM curves: basin x energy tertile ---
pt[, E_tertile := cut(energy,
  breaks = quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  labels = c("low", "mid", "high"), include.lowest = TRUE)]

km_panels <- list()
for (b in basins) {
  sub <- pt[basin_id == b & !is.na(E_tertile)]
  if (nrow(sub) < 30) next

  sf <- survfit(Surv(time_years, death_bin) ~ E_tertile, data = sub)

  # Extract KM data manually (no survminer dependency)
  km_data <- data.table(
    time = sf$time,
    surv = sf$surv,
    strata = rep(names(sf$strata), sf$strata)
  )
  km_data[, E_tertile := sub("E_tertile=", "", km_data$strata)]
  km_data[, basin_label := sprintf("Basin %d (n=%s)",
                                    b, formatC(nrow(sub), big.mark = ","))]

  # Add HR annotation
  hr_row <- cox_dt[basin_id == b]
  if (!is.na(hr_row$HR)) {
    km_data[, hr_label := sprintf("HR=%.3f [%.3f,%.3f] p=%.1e",
                                   hr_row$HR, hr_row$HR_lower,
                                   hr_row$HR_upper, hr_row$pvalue)]
  } else {
    km_data[, hr_label := "HR: N/A"]
  }

  km_panels[[length(km_panels) + 1]] <- km_data
}

if (length(km_panels) > 0) {
  km_all <- rbindlist(km_panels)

  # Add origin point for each group
  origins <- unique(km_all[, .(basin_label, E_tertile, hr_label)])
  origins[, `:=`(time = 0, surv = 1)]
  km_all <- rbind(origins, km_all, fill = TRUE)

  g <- ggplot(km_all, aes(x = time, y = surv, color = E_tertile)) +
    geom_step(linewidth = 0.8) +
    facet_wrap(~ basin_label, scales = "free_x") +
    scale_color_manual(
      values = c(low = "#d94a3a", mid = "#f0a830", high = "#3a6dd9"),
      labels = c(low = "Low energy (severe)",
                 mid = "Mid energy",
                 high = "High energy (mild)")
    ) +
    geom_text(data = unique(km_all[, .(basin_label, hr_label)]),
              aes(label = hr_label),
              x = Inf, y = 0.02, hjust = 1.05, vjust = 0,
              size = 3, color = "gray30", inherit.aes = FALSE) +
    labs(title = sprintf("%s: KM survival by energy tertile within each basin",
                         disease),
         x = "Years", y = "Survival probability", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  w <- max(8, 5 * n_basins)
  h <- max(5, 4 * ceiling(n_basins / 3))
  ggsave(file.path(out_dir, sprintf("km_basin_energy_%s.png", disease)),
         g, width = min(w, 20), height = min(h, 16), dpi = 150)
  cat(sprintf("  -> km_basin_energy_%s.png\n", disease))
}

cat("Done\n")
