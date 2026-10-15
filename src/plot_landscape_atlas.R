## Landscape Atlas Visualization
##
## Generates per-setting landscape summary figures and a cross-setting
## atlas heatmap for comparing energy direction across disease families
## and preprocessing strategies.
##
## Usage (per-setting):
##   Rscript src/plot_landscape_atlas.R \
##     <person_tsv> <cox_tsv> <basin_tsv> <energy_tsv> \
##     <output_dir> <disease> <preproc> <dim>
##
## Usage (atlas heatmap):
##   Rscript src/plot_landscape_atlas.R \
##     --atlas-heatmap <summary_tsv> <output_png>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

# -------------------------------------------------------------------
# Atlas heatmap mode
# -------------------------------------------------------------------
if (length(args) >= 1 && args[1] == "--atlas-heatmap") {
  summary_path <- args[2]
  out_path <- args[3]

  cat("=== Atlas Heatmap ===\n")
  sm <- fread(summary_path)

  sm[, direction_label := ifelse(
    is.na(direction_ok), "N/A",
    ifelse(direction_ok == TRUE, "PASS", "FAIL")
  )]
  sm[, hr_label := sprintf("%.3f", HR)]

  sm[, disease_family := factor(disease_family,
       levels = rev(sort(unique(disease_family))))]

  g <- ggplot(sm, aes(x = outcome, y = disease_family,
                       fill = direction_label)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = hr_label), size = 3) +
    scale_fill_manual(
      values = c(PASS = "#4CAF50", FAIL = "#F44336", "N/A" = "gray80"),
      name = "Energy direction"
    ) +
    labs(title = "Energy Cox Sanity Check Atlas",
         subtitle = "HR per +1 SD Energy. PASS = HR < 1 (low energy = high risk)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(out_path, g, width = 8, height = 16, dpi = 200)
  cat(sprintf("Written: %s\n", out_path))
  quit(save = "no")
}

# -------------------------------------------------------------------
# Per-setting mode
# -------------------------------------------------------------------
stopifnot(length(args) == 8)

person_path <- args[1]
cox_path    <- args[2]
basin_path  <- args[3]
energy_path <- args[4]
out_dir     <- args[5]
disease     <- args[6]
preproc     <- args[7]
dim_val     <- args[8]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== plot_landscape_atlas: %s / %s / dim=%s ===\n",
            disease, preproc, dim_val))

# --- Load data ---
pt <- fread(person_path)
cox_res <- fread(cox_path)
energy_data <- fread(energy_path, header = FALSE)
setnames(energy_data, "energy")
energy_data[, pattern_id := .I]

# --- Summary figure: energy distribution + basin assignment ---
if ("energy" %in% names(pt) && "basin_id" %in% names(pt)) {
  pt[, basin_label := paste0("Basin ", basin_id)]

  g1 <- ggplot(pt[is.finite(energy)], aes(x = energy, fill = basin_label)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(title = sprintf("%s / %s / dim=%s", disease, preproc, dim_val),
         subtitle = sprintf("n=%d patients", nrow(pt)),
         x = "Energy", y = "Count", fill = NULL) +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "landscape_summary.png"), g1,
         width = 8, height = 5, dpi = 200)

  # --- KM by energy tertile ---
  time_col <- intersect(c("time_years", "obs_years"), names(pt))[1]
  event_col <- intersect(c("event_death", "death"), names(pt))[1]
  if (!is.na(time_col) && !is.na(event_col)) {
    library(survival)

    if (is.character(pt[[event_col]])) {
      pt[, (event_col) := as.integer(get(event_col) == "dead")]
    }

    pt_valid <- pt[is.finite(get(time_col)) & get(time_col) > 0 & is.finite(energy)]
    if (nrow(pt_valid) > 0 && sum(pt_valid[[event_col]], na.rm = TRUE) >= 5) {
      pt_valid[, E_tertile := cut(energy,
        breaks = quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
        labels = c("low", "mid", "high"),
        include.lowest = TRUE)]

      km <- survfit(
        Surv(pt_valid[[time_col]], pt_valid[[event_col]]) ~ E_tertile,
        data = pt_valid)
      km_df <- data.frame(
        time = km$time, surv = km$surv,
        lower = km$lower, upper = km$upper,
        E_tertile = rep(names(km$strata), km$strata)
      )
      km_df$E_tertile <- factor(
        sub("E_tertile=", "", km_df$E_tertile),
        levels = c("low", "mid", "high"))

      g2 <- ggplot(km_df, aes(x = time, y = surv,
                               color = E_tertile, fill = E_tertile)) +
        geom_step(linewidth = 1.0) +
        geom_ribbon(aes(ymin = lower, ymax = upper),
                    alpha = 0.15, color = NA) +
        scale_color_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                       high = "#d94a3a")) +
        scale_fill_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                      high = "#d94a3a")) +
        labs(title = sprintf("KM by energy tertile: %s / %s",
                             disease, preproc),
             x = "Years", y = "Survival probability",
             color = "Energy", fill = "Energy") +
        theme_minimal(base_size = 14)
      ggsave(file.path(out_dir, "km_energy.png"), g2,
             width = 8, height = 5, dpi = 200)
    }
  }
}

# --- Forest plot of energy HR ---
if (nrow(cox_res) > 0 && any(!is.na(cox_res$HR))) {
  cox_res[, outcome := factor(outcome, levels = c("mortality", "hospitalization"))]
  cox_res[, direction := ifelse(direction_ok, "PASS", "FAIL")]

  g3 <- ggplot(cox_res[!is.na(HR)],
               aes(x = HR, y = outcome, color = direction)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray40") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_point(size = 4) +
    scale_color_manual(values = c(PASS = "#4CAF50", FAIL = "#F44336")) +
    labs(title = sprintf("Energy HR sanity check: %s / %s / dim=%s",
                         disease, preproc, dim_val),
         x = "HR per +1 SD Energy (95% CI)", y = NULL, color = NULL) +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "forest_energy.png"), g3,
         width = 8, height = 4, dpi = 200)
}

cat(sprintf("Done. Figures in %s\n", out_dir))
