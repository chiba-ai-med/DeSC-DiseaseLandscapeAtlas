## KM survival curves per basin
##
## Produces two plots:
##   1. km_inter_basin.png — Basin comparison (Surv ~ basin_id)
##   2. km_per_basin.png   — Within-basin energy tertile (faceted)
##
## Usage:
##   Rscript src/km_per_basin.R <disease> <patient_scores_tsv> \
##     <tekiyo_csv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
})

args <- commandArgs(trailingOnly = TRUE)
disease     <- args[1]
scores_path <- args[2]
tekiyo_path <- args[3]
out_dir     <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

disease_label <- gsub("_", " ", paste0(toupper(substring(disease, 1, 1)),
                                        substring(disease, 2)))
cat(sprintf("=== KM per basin: %s ===\n", disease))

pt <- fread(scores_path)

tek <- fread(tekiyo_path,
             select = c("kojin_id", "birth_ym", "sex_code",
                         "observable_start_ym", "observable_end_ym",
                         "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
pt[, time_years := parse_ym(observable_end_ym) - parse_ym(observable_start_ym)]
pt[, event_death := as.integer(!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1)]

pt <- pt[is.finite(time_years) & time_years > 0 & is.finite(energy)]

basins <- sort(unique(pt[basin_id > 0]$basin_id))
n_basins <- length(basins)
cat(sprintf("  %d patients, %d basins\n", nrow(pt), n_basins))

## ── 1. Inter-basin KM (Basin membership comparison) ──
if (n_basins >= 2) {
  pt_ib <- pt[basin_id > 0]
  pt_ib[, Basin := sprintf("Basin %d (n=%s)", basin_id,
                            format(.N, big.mark = ",")), by = basin_id]
  basin_order <- pt_ib[, .(bid = basin_id[1]), by = Basin][order(bid)]$Basin
  pt_ib[, Basin := factor(Basin, levels = basin_order)]

  fit_ib <- survfit(Surv(time_years, event_death) ~ Basin, data = pt_ib)

  n_colors <- n_basins
  basin_colors <- RColorBrewer::brewer.pal(max(3, n_colors), "Set1")[seq_len(n_colors)]
  basin_labels <- levels(pt_ib$Basin)

  g_ib <- ggsurvplot(fit_ib, data = pt_ib,
    pval = TRUE, pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.25,
    tables.y.text = TRUE,
    palette = basin_colors,
    xlab = "Years", ylab = "Survival probability",
    title = sprintf("%s - Survival by basin", disease_label),
    legend.title = "",
    legend.labs = basin_labels,
    ggtheme = theme_bw(base_size = 14),
    tables.theme = theme_bw(base_size = 12),
    font.main = c(16, "bold"),
    font.x = 13, font.y = 13,
    font.tickslab = 11,
    risk.table.fontsize = 3.5,
    surv.plot.height = 0.72
  )

  g_ib$plot <- g_ib$plot + theme(
    axis.title.y = element_text(margin = margin(r = 0)),
    plot.margin = margin(5, 15, 0, 2))
  g_ib$table <- g_ib$table + theme(
    plot.margin = margin(0, 15, 5, 2))

  out_ib <- file.path(out_dir, "km_inter_basin.png")
  png(out_ib, width = 10, height = 7, units = "in", res = 300, bg = "white")
  print(g_ib)
  dev.off()
  cat(sprintf("  Written: %s\n", out_ib))
}

## ── 2. Per-basin KM (energy tertile within each basin) ──
plot_list <- list()
for (b in basins) {
  sub <- pt[basin_id == b]
  if (nrow(sub) < 30 || sum(sub$event_death) < 10) next

  brks <- quantile(sub$energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (length(unique(brks)) < 4) next

  sub[, E_tertile := cut(energy, breaks = brks,
                         labels = c("Low", "Mid", "High"),
                         include.lowest = TRUE)]
  sub <- sub[!is.na(E_tertile)]
  if (nrow(sub) < 30) next

  fit <- survfit(Surv(time_years, event_death) ~ E_tertile, data = sub)

  g <- ggsurvplot(fit, data = sub,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.30,
    palette = c("#3a6dd9", "gray50", "#d94a3a"),
    xlab = "Years", ylab = "Survival probability",
    title = sprintf("%s - Basin %d (n=%s)",
                    disease_label, b, format(nrow(sub), big.mark = ",")),
    legend.title = "Energy",
    ggtheme = theme_bw(base_size = 13),
    font.main = c(15, "bold"),
    font.x = 12, font.y = 12,
    font.tickslab = 10,
    risk.table.fontsize = 3.2
  )

  out_path <- file.path(out_dir, sprintf("km_energy_basin%d.png", b))
  png(out_path, width = 8, height = 6, units = "in", res = 300, bg = "white")
  print(g)
  dev.off()
  cat(sprintf("  Written: %s\n", out_path))
  plot_list[[length(plot_list) + 1]] <- g
}

## Combined faceted version
if (length(plot_list) >= 2) {
  n_col <- min(2, length(plot_list))
  n_row <- ceiling(length(plot_list) / n_col)
  g_arr <- arrange_ggsurvplots(plot_list, ncol = n_col, nrow = n_row)

  out_combined <- file.path(out_dir, "km_per_basin.png")
  png(out_combined, width = 8 * n_col, height = 6 * n_row,
      units = "in", res = 200, bg = "white")
  print(g_arr)
  dev.off()
  cat(sprintf("  Written: %s\n", out_combined))
}

cat("Done.\n")
