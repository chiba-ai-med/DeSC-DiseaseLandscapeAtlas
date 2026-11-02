## Generate landmark-design KM for diabetes (Fig 6b candidate)
##
## Reads landmark_cohort.tsv (fu_time, event_death) instead of
## patient_scores.tsv (observable_end_ym - observable_start_ym).
## Style matches fig_casestudy_panels.R Panel b.
##
## Usage: Rscript paper/claude/review_candidates/fig06b_landmark_km.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

cohort_path <- "output/diabetes/all_comorbid/landmark/landmark_cohort.tsv"
out_dir <- "paper/claude/review_candidates"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Landmark KM for diabetes (Fig 6b candidate) ===\n")

# ── Load and filter ──
dt <- fread(cohort_path)
cat(sprintf("  Total rows: %d\n", nrow(dt)))

# Mortality-eligible: death flag available, valid basin, positive follow-up
dt <- dt[shibou_flg_riyouka == 1L & basin_id > 0L & fu_time > 0]
cat(sprintf("  After filter (riyouka=1, basin>0, fu>0): %d\n", nrow(dt)))

# ── Basin summary ──
cat("\n  Basin summary:\n")
for (b in sort(unique(dt$basin_id))) {
  sub <- dt[basin_id == b]
  cat(sprintf("    Basin %d: n=%s, deaths=%s (%.1f%%), median fu=%.2f yr\n",
              b, format(nrow(sub), big.mark = ","),
              format(sum(sub$event_death), big.mark = ","),
              100 * mean(sub$event_death),
              median(sub$fu_time)))
}

# ── Prepare for KM ──
dt[, Basin := sprintf("Basin %d", basin_id)]
basin_order <- sort(unique(dt$Basin))
dt[, Basin := factor(Basin, levels = basin_order)]

n_basins <- length(basin_order)
basin_colors <- brewer.pal(max(3, n_basins), "Set1")[seq_len(n_basins)]
names(basin_colors) <- basin_order

# ── Fit KM ──
fit_km <- survfit(Surv(fu_time, event_death) ~ Basin, data = dt)

survfit_to_dt <- function(fit) {
  strata_names <- names(fit$strata)
  strata_lens <- fit$strata
  data.table(time = fit$time, surv = fit$surv,
             upper = fit$upper, lower = fit$lower,
             n.risk = fit$n.risk,
             strata = rep(strata_names, strata_lens))
}

sdf <- survfit_to_dt(fit_km)
sdf[, strata := gsub("Basin=", "", strata)]
time0 <- data.table(time = 0, surv = 1, upper = 1, lower = 1,
                    n.risk = NA_real_, strata = unique(sdf$strata))
sdf <- rbind(time0, sdf)

# ── Log-rank test ──
sd <- survdiff(Surv(fu_time, event_death) ~ Basin, data = dt)
pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
pv_label <- if (pv < 0.0001) "Log-rank\np < 0.0001" else sprintf("Log-rank\np = %.4f", pv)

# ── Plot (matching fig_casestudy_panels.R style) ──
ksz <- 14
g_km <- ggplot(sdf, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12,
              stat = "identity", outline.type = "upper") +
  scale_color_manual(values = basin_colors, name = "") +
  scale_fill_manual(values = basin_colors, name = "") +
  annotate("text", x = 2.5, y = 0.15, label = pv_label,
           size = ksz * 3 / .pt, fontface = "bold", hjust = 0) +
  labs(title = "Diabetes - Survival by basin (landmark)",
       x = "Years from landmark", y = "Survival probability") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 8)) +
  theme_bw(base_size = ksz) +
  theme(plot.title = element_text(size = ksz * 3, face = "bold"),
        axis.text = element_text(size = ksz * 3),
        axis.title = element_text(size = ksz * 3),
        legend.text = element_text(size = ksz * 2.5),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        plot.background = element_rect(fill = "white", colour = NA))

# ── Risk table ──
risk_times <- seq(0, 8, by = 2)
risk_dt <- data.table()
for (b in basin_order) {
  sub <- dt[Basin == b]
  for (tt in risk_times) {
    risk_dt <- rbind(risk_dt, data.table(
      time = tt, Basin = b, n_risk = sum(sub$fu_time >= tt)))
  }
}
risk_dt[, Basin := factor(Basin, levels = rev(basin_order))]
risk_dt[, label := format(n_risk, big.mark = ",")]

g_risk <- ggplot(risk_dt, aes(x = time, y = Basin, label = label, color = Basin)) +
  geom_text(size = ksz * 2.5 / .pt, fontface = "bold") +
  scale_color_manual(values = basin_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.5, 8), breaks = risk_times) +
  labs(x = NULL, y = NULL, title = "Number at risk") +
  theme_bw(base_size = ksz) +
  theme(plot.title = element_text(size = ksz * 2, face = "bold"),
        axis.text.x = element_text(size = ksz * 3),
        axis.text.y = element_text(size = ksz * 2),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 10, 10, 80))

# ── Compose and save ──
out_path <- file.path(out_dir, "Fig06_b_diabetes_km_landmark.png")
png(out_path, width = 24, height = 18, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.75, 0.25), "npc"))))
print(g_km, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g_risk, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

cat(sprintf("\n  Written: %s\n", out_path))
cat("  Done.\n")
