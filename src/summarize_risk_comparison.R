## Cross-disease summary of risk score comparison
##
## Reads per-disease risk_score_comparison.tsv files and produces:
## - Summary table
## - Heatmap of C-statistics
## - Scatter: Energy C vs CCI C
##
## Usage:
##   Rscript src/summarize_risk_comparison.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cat("=== Risk score comparison summary ===\n")

out_dir <- "summary/risk_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Collect all per-disease results
files <- Sys.glob("output/*/all_comorbid/comparison/risk_score_comparison.tsv")
cat(sprintf("  Found %d disease files\n", length(files)))

all_dt <- rbindlist(lapply(files, fread))

# --- Summary table ---
fwrite(all_dt, file.path(out_dir, "all_comparisons.tsv"), sep = "\t")

# --- For each endpoint, compare models ---
for (ep in unique(all_dt$endpoint)) {
  dt_ep <- all_dt[endpoint == ep]

  # Pivot: disease × model -> C
  wide <- dcast(dt_ep, disease ~ model, value.var = "C")
  fwrite(wide, file.path(out_dir, sprintf("c_statistics_%s.tsv", ep)), sep = "\t")

  # How often does Energy beat CCI / Elixhauser?
  if ("Energy" %in% names(wide) && "CCI" %in% names(wide)) {
    n_energy_beats_cci <- sum(wide$Energy > wide$CCI, na.rm = TRUE)
    n_total <- sum(!is.na(wide$Energy) & !is.na(wide$CCI))
    cat(sprintf("  [%s] Energy > CCI: %d / %d diseases (%.0f%%)\n",
        ep, n_energy_beats_cci, n_total, 100 * n_energy_beats_cci / n_total))
  }
  if ("Energy" %in% names(wide) && "Elixhauser" %in% names(wide)) {
    n_energy_beats_elix <- sum(wide$Energy > wide$Elixhauser, na.rm = TRUE)
    n_total <- sum(!is.na(wide$Energy) & !is.na(wide$Elixhauser))
    cat(sprintf("  [%s] Energy > Elixhauser: %d / %d diseases (%.0f%%)\n",
        ep, n_energy_beats_elix, n_total, 100 * n_energy_beats_elix / n_total))
  }

  # --- Heatmap ---
  core_models <- c("Energy", "CCI", "Elixhauser", "N_codes",
                    "Energy+Age+Sex", "CCI+Age+Sex", "Elix+Age+Sex",
                    "Basin+Age+Sex", "Energy+CCI+Age+Sex", "Energy+Elix+Age+Sex",
                    "Basin+CCI+Age+Sex", "Energy+Basin+CCI+Age+Sex")
  dt_plot <- dt_ep[model %in% core_models]
  dt_plot[, model := factor(model, levels = core_models)]

  g <- ggplot(dt_plot, aes(x = model, y = disease, fill = C)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", C)), size = 2.2) +
    scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                         midpoint = 0.5, limits = c(0.4, 0.8),
                         oob = scales::squish, name = "C-statistic") +
    labs(title = sprintf("Risk score comparison: %s", ep),
         x = "", y = "") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))

  ggsave(file.path(out_dir, sprintf("heatmap_%s.png", ep)),
         g, width = 14, height = 16, dpi = 200, bg = "white")

  # --- Scatter: Energy vs CCI C-statistic ---
  if ("Energy" %in% names(wide) && "CCI" %in% names(wide)) {
    g2 <- ggplot(wide, aes(x = CCI, y = Energy, label = disease)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(color = "steelblue", size = 2.5) +
      geom_text(size = 2.5, hjust = -0.1, vjust = -0.3) +
      labs(title = sprintf("C-statistic: Energy vs CCI (%s)", ep),
           subtitle = "Points above diagonal = Energy is better",
           x = "C-statistic (CCI)", y = "C-statistic (Energy)") +
      coord_equal() +
      theme_minimal(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(out_dir, sprintf("scatter_energy_vs_cci_%s.png", ep)),
           g2, width = 10, height = 10, dpi = 200, bg = "white")
  }

  # --- Scatter: Energy vs Elixhauser ---
  if ("Energy" %in% names(wide) && "Elixhauser" %in% names(wide)) {
    g3 <- ggplot(wide, aes(x = Elixhauser, y = Energy, label = disease)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(color = "coral", size = 2.5) +
      geom_text(size = 2.5, hjust = -0.1, vjust = -0.3) +
      labs(title = sprintf("C-statistic: Energy vs Elixhauser (%s)", ep),
           subtitle = "Points above diagonal = Energy is better",
           x = "C-statistic (Elixhauser)", y = "C-statistic (Energy)") +
      coord_equal() +
      theme_minimal(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(out_dir, sprintf("scatter_energy_vs_elix_%s.png", ep)),
           g3, width = 10, height = 10, dpi = 200, bg = "white")
  }

  # --- Bar: Energy+CCI vs CCI alone (added value) ---
  if ("Energy+CCI+Age+Sex" %in% names(wide) && "CCI+Age+Sex" %in% names(wide)) {
    wide[, delta_cci := `Energy+CCI+Age+Sex` - `CCI+Age+Sex`]
    setorder(wide, -delta_cci)
    wide[, disease_f := factor(disease, levels = disease)]
    g4 <- ggplot(wide, aes(x = disease_f, y = delta_cci)) +
      geom_col(aes(fill = delta_cci > 0), show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "coral")) +
      coord_flip() +
      labs(title = sprintf("Added value of Energy over CCI+Age+Sex (%s)", ep),
           subtitle = "Delta C-statistic: (Energy+CCI+Age+Sex) - (CCI+Age+Sex)",
           x = "", y = "Delta C-statistic") +
      theme_minimal(base_size = 10) +
      theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(out_dir, sprintf("bar_delta_cci_%s.png", ep)),
           g4, width = 10, height = 14, dpi = 200, bg = "white")
  }
}

cat("Done.\n")
