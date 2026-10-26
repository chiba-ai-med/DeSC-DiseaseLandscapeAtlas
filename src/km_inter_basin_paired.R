## Generate paired basin-comparison KM plots (mortality + hospitalization)
##
## For a disease with ≥2 basins, produces a side-by-side figure:
##   (a) Mortality: Surv(obs_years, death) ~ Basin
##   (b) Hospitalization: Surv(time_to_first_admission, event_hosp) ~ Basin
##
## Usage:
##   Rscript src/km_inter_basin_paired.R <disease> <base_dir> <hosp_specific> <tekiyo> <out_png>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
disease   <- args[1]
base_dir  <- args[2]
hosp_path <- args[3]
tek_path  <- args[4]
out_path  <- args[5]

cat(sprintf("=== KM inter-basin paired: %s ===\n", disease))

pt <- fread(file.path(base_dir, "patient/patient_scores.tsv"))
pt[, kojin_id := as.character(kojin_id)]

basins <- sort(unique(pt[basin_id > 0]$basin_id))
n_basins <- length(basins)
if (n_basins < 2) {
  cat("Fewer than 2 basins, skipping.\n")
  quit(save = "no")
}

tek <- fread(tek_path,
             select = c("kojin_id", "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
pt[, time_years := parse_ym(observable_end_ym) - parse_ym(observable_start_ym)]
pt[, event_death := as.integer(!is.na(shibou_flg) & shibou_flg == 1)]

cat(sprintf("Before riyouka filter: %d patients\n", nrow(pt)))
pt <- pt[!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1]
cat(sprintf("After riyouka filter: %d patients\n", nrow(pt)))

pt_death <- pt[basin_id > 0 & is.finite(time_years) & time_years > 0]
pt_death[, Basin := factor(sprintf("Basin %d", basin_id))]

n_colors <- n_basins
basin_colors <- brewer.pal(max(3, n_colors), "Set1")[seq_len(n_colors)]

survfit_to_df <- function(fit) {
  strata_names <- names(fit$strata)
  strata_lens <- fit$strata
  strata_vec <- rep(strata_names, strata_lens)
  data.table(
    time = fit$time,
    surv = fit$surv,
    upper = fit$upper,
    lower = fit$lower,
    n.risk = fit$n.risk,
    strata = strata_vec
  )
}

make_km_plot <- function(fit, df, title, ylab, colors, event_col) {
  sdf <- survfit_to_df(fit)
  sdf[, strata := gsub("Basin=", "", strata)]

  time0 <- data.table(
    time = 0, surv = 1, upper = 1, lower = 1, n.risk = NA_real_,
    strata = unique(sdf$strata)
  )
  sdf <- rbind(time0, sdf)

  p <- ggplot(sdf, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_step(linewidth = 0.8) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, stat = "identity",
                outline.type = "upper") +
    scale_color_manual(values = setNames(colors, unique(sdf$strata))) +
    scale_fill_manual(values = setNames(colors, unique(sdf$strata))) +
    labs(x = "Years", y = ylab, title = title, color = "", fill = "") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    ) +
    coord_cartesian(ylim = c(0, 1))

  basin_labels <- unique(sdf$strata)
  n_per_basin <- sapply(basin_labels, function(b) {
    nrow(df[Basin == b])
  })
  event_per_basin <- sapply(basin_labels, function(b) {
    sum(df[Basin == b][[event_col]])
  })
  rate_per_basin <- round(event_per_basin / n_per_basin * 100, 1)
  legend_text <- paste0(basin_labels, " (n=", formatC(n_per_basin, format = "d", big.mark = ","),
                        "; ", rate_per_basin, "%)")
  p <- p + scale_color_manual(values = setNames(colors, basin_labels),
                              labels = setNames(legend_text, basin_labels)) +
           scale_fill_manual(values = setNames(colors, basin_labels),
                             labels = setNames(legend_text, basin_labels))
  p
}

fit_death <- survfit(Surv(time_years, event_death) ~ Basin, data = pt_death)
g_death <- make_km_plot(fit_death, pt_death, "(a) Mortality",
                        "Survival probability", basin_colors, "event_death")

hosp <- fread(hosp_path)
hosp[, kojin_id := as.character(kojin_id)]
pt_hosp <- merge(pt[basin_id > 0], hosp[, .(kojin_id, first_admission_ymd)],
                 by = "kojin_id", all.x = TRUE)

parse_ymd <- function(s) {
  y <- as.numeric(substr(s, 1, 4))
  m <- as.numeric(substr(s, 6, 7))
  d <- ifelse(nchar(s) >= 10, as.numeric(substr(s, 9, 10)), 15)
  y + (m - 1 + (d - 1) / 30) / 12
}

pt_hosp[, start_y := parse_ym(observable_start_ym)]
pt_hosp[, end_y := parse_ym(observable_end_ym)]
pt_hosp[, adm_y := parse_ymd(first_admission_ymd)]
pt_hosp[, has_adm := !is.na(adm_y) & adm_y >= start_y & adm_y <= end_y]
pt_hosp[, time_hosp := ifelse(has_adm, adm_y - start_y, end_y - start_y)]
pt_hosp[, event_hosp := as.integer(has_adm)]

pt_hosp <- pt_hosp[is.finite(time_hosp) & time_hosp > 0]
pt_hosp[, Basin := factor(sprintf("Basin %d", basin_id))]

fit_hosp <- survfit(Surv(time_hosp, event_hosp) ~ Basin, data = pt_hosp)
g_hosp <- make_km_plot(fit_hosp, pt_hosp, "(b) Hospitalization",
                       "Hospitalization-free probability", basin_colors, "event_hosp")

disease_label <- gsub("_", " ", disease)
disease_label <- paste0(toupper(substr(disease_label, 1, 1)),
                        substr(disease_label, 2, nchar(disease_label)))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
png(out_path, width = 16, height = 7, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.06, 0.94), "npc"))))
grid.text(disease_label, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          gp = gpar(fontsize = 18, fontface = "bold"))
print(g_death, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(g_hosp, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

cat(sprintf("Written: %s\n", out_path))
cat(sprintf("  Mortality: %s\n", paste(sprintf("Basin %d (n=%d)",
    basins, table(pt_death$Basin)), collapse = ", ")))
cat(sprintf("  Hospitalization: %s\n", paste(sprintf("Basin %d (n=%d)",
    basins, table(pt_hosp$Basin)), collapse = ", ")))
cat("Done.\n")
