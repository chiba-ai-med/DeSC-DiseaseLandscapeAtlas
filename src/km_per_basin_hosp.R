## KM curves per basin, stratified by energy tertile — hospitalization endpoint
##
## Same structure as km_per_basin.R but event = first admission.
##
## Usage:
##   Rscript src/km_per_basin_hosp.R <disease> <patient_scores_tsv> \
##     <hosp_tsv> <tekiyo_csv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease     <- args[1]
scores_path <- args[2]
hosp_path   <- args[3]
tekiyo_path <- args[4]
out_dir     <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== KM per basin (hospitalization): %s ===\n", disease))

pt <- fread(scores_path)

tek <- fread(tekiyo_path,
             select = c("kojin_id", "birth_ym", "sex_code",
                         "observable_start_ym", "observable_end_ym"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

hosp <- fread(hosp_path)
hosp[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, hosp[, .(kojin_id, first_admission_ymd)],
            by = "kojin_id", all.x = TRUE)

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
parse_ymd <- function(s) {
  y <- as.numeric(substr(s, 1, 4))
  m <- as.numeric(substr(s, 6, 7))
  d <- ifelse(nchar(s) >= 10, as.numeric(substr(s, 9, 10)), 15)
  y + (m - 1 + (d - 1) / 30) / 12
}

pt[, start_y := parse_ym(observable_start_ym)]
pt[, end_y := parse_ym(observable_end_ym)]
pt[, adm_y := parse_ymd(first_admission_ymd)]
pt[, has_adm := !is.na(adm_y) & adm_y >= start_y & adm_y <= end_y]
pt[, time_years := ifelse(has_adm, adm_y - start_y, end_y - start_y)]
pt[, event_hosp := as.integer(has_adm)]

pt <- pt[is.finite(time_years) & time_years > 0 & is.finite(energy)]

basins <- sort(unique(pt$basin_id))
n_basins <- length(basins)
cat(sprintf("  %d patients, %d basins, %d hospitalized\n",
            nrow(pt), n_basins, sum(pt$event_hosp)))

# --- Per-basin KM (energy tertile within each basin) ---
km_list <- list()
for (b in basins) {
  sub <- pt[basin_id == b]
  if (nrow(sub) < 30 || sum(sub$event_hosp) < 10) next

  brks <- quantile(sub$energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (length(unique(brks)) < 4) next

  sub[, E_tertile := cut(energy, breaks = brks,
                         labels = c("low", "mid", "high"),
                         include.lowest = TRUE)]
  sub <- sub[!is.na(E_tertile)]
  if (nrow(sub) < 30) next

  km <- survfit(Surv(time_years, event_hosp) ~ E_tertile, data = sub)
  km_df <- data.frame(
    time = km$time, surv = km$surv,
    lower = km$lower, upper = km$upper,
    strata = rep(names(km$strata), km$strata)
  )
  km_df$E_tertile <- factor(
    sub("E_tertile=", "", km_df$strata),
    levels = c("low", "mid", "high"))
  km_df$basin <- paste0(disease, " / Basin ", b,
                        " (n=", format(nrow(sub), big.mark=","), ")")
  km_list[[length(km_list) + 1]] <- km_df
}

out_path <- file.path(out_dir, "km_per_basin_hosp.png")
if (length(km_list) > 0) {
  km_all <- do.call(rbind, km_list)

  height <- 4 * ceiling(n_basins / 2)
  width <- min(14, 7 * min(n_basins, 2))

  g <- ggplot(km_all, aes(x = time, y = surv,
                           color = E_tertile, fill = E_tertile)) +
    geom_step(linewidth = 0.8) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.12, color = NA) +
    facet_wrap(~ basin, scales = "free_x") +
    scale_color_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                   high = "#d94a3a")) +
    scale_fill_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                  high = "#d94a3a")) +
    labs(x = "Years", y = "Hospitalization-free probability",
         color = "Energy", fill = "Energy") +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(face = "bold"))

  ggsave(out_path, g, width = width, height = height, dpi = 200)
  cat(sprintf("  Written: %s\n", out_path))
} else {
  # No basin had enough variation for tertile split; write placeholder
  png(out_path, width = 400, height = 200)
  plot.new(); text(0.5, 0.5, paste(disease, "\nNo energy tertile variation"), cex = 1.2)
  dev.off()
  cat(sprintf("  Written (placeholder): %s\n", out_path))
}

# --- Combined: all basins on one plot ---
if (n_basins >= 2) {
  pt_valid <- pt[is.finite(energy) & is.finite(time_years) & time_years > 0]
  pt_valid[, E_tertile := {
    brks <- quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    if (length(unique(brks)) < 4) NA_character_
    else as.character(cut(energy, breaks = brks,
      labels = c("low", "mid", "high"), include.lowest = TRUE))
  }, by = basin_id]
  pt_valid <- pt_valid[!is.na(E_tertile)]
  pt_valid[, E_tertile := factor(E_tertile, levels = c("low", "mid", "high"))]
  pt_valid[, group := paste0("Basin", basin_id, "_", E_tertile)]

  if (nrow(pt_valid) < 30 || length(unique(pt_valid$group)) < 2) {
    cat("  Combined plot: too few data after tertile split, skipping.\n")
  } else {
  km2 <- survfit(Surv(time_years, event_hosp) ~ group, data = pt_valid)
  km2_df <- data.frame(
    time = km2$time, surv = km2$surv,
    lower = km2$lower, upper = km2$upper,
    strata = rep(names(km2$strata), km2$strata)
  )
  km2_df$group <- sub("group=", "", km2_df$strata)

  g2 <- ggplot(km2_df, aes(x = time, y = surv,
                            color = group, linetype = group)) +
    geom_step(linewidth = 0.7) +
    labs(x = "Years", y = "Hospitalization-free probability",
         color = NULL, linetype = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")

  out_path2 <- file.path(out_dir, "km_basin_energy_combined_hosp.png")
  ggsave(out_path2, g2, width = 10, height = 6, dpi = 200)
  cat(sprintf("  Written: %s\n", out_path2))
  }
}

cat("Done.\n")
