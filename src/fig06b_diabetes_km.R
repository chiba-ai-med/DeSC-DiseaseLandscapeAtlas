suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"

scores_path <- file.path(base, "output/diabetes/all_comorbid/patient/patient_scores.tsv")
tekiyo_path <- file.path(base, "data/tekiyo.csv")

pt <- fread(scores_path)
tek <- fread(tekiyo_path,
             select = c("kojin_id", "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
pt[, time_years := parse_ym(observable_end_ym) - parse_ym(observable_start_ym)]
pt[, event_death := as.integer(!is.na(shibou_flg) & shibou_flg == 1)]

pt <- pt[!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1]
pt <- pt[is.finite(time_years) & time_years > 0 & is.finite(energy)]

pt_ib <- pt[basin_id > 0]
pt_ib[, Basin := sprintf("Basin %d", basin_id)]
basin_order <- pt_ib[, .(bid = basin_id[1]), by = Basin][order(bid)]$Basin
pt_ib[, Basin := factor(Basin, levels = basin_order)]

fit_ib <- survfit(Surv(time_years, event_death) ~ Basin, data = pt_ib)

n_basins <- length(unique(pt_ib$basin_id))
basin_colors <- brewer.pal(max(3, n_basins), "Set1")[seq_len(n_basins)]
names(basin_colors) <- basin_order

survfit_to_dt <- function(fit) {
  strata_names <- names(fit$strata)
  strata_lens <- fit$strata
  strata_vec <- rep(strata_names, strata_lens)
  data.table(
    time = fit$time, surv = fit$surv,
    upper = fit$upper, lower = fit$lower,
    n.risk = fit$n.risk, strata = strata_vec
  )
}

sdf <- survfit_to_dt(fit_ib)
sdf[, strata := gsub("Basin=", "", strata)]
time0 <- data.table(time = 0, surv = 1, upper = 1, lower = 1,
                    n.risk = NA_real_, strata = unique(sdf$strata))
sdf <- rbind(time0, sdf)

# Log-rank test
sd <- survdiff(Surv(time_years, event_death) ~ Basin, data = pt_ib)
pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
pv_label <- if (pv < 0.0001) "Log-rank\np < 0.0001" else sprintf("Log-rank\np = %.4f", pv)

sz <- 14  # base, then 3x for axis/logrank

g_km <- ggplot(sdf, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12,
              stat = "identity", outline.type = "upper") +
  scale_color_manual(values = basin_colors, name = "") +
  scale_fill_manual(values = basin_colors, name = "") +
  annotate("text", x = 2.5, y = 0.15, label = pv_label,
           size = sz * 3 / .pt, fontface = "bold", hjust = 0) +
  labs(title = "Diabetes - Survival by basin",
       x = "Years", y = "Survival probability") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 8)) +
  theme_bw(base_size = sz) +
  theme(
    plot.title = element_text(size = sz * 3, face = "bold"),
    axis.text = element_text(size = sz * 3),
    axis.title = element_text(size = sz * 3),
    legend.text = element_text(size = sz * 2.5),
    legend.position = "top",
    legend.key.width = unit(2, "cm"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Risk table
risk_times <- seq(0, 8, by = 2)
risk_dt <- data.table()
for (b in basin_order) {
  sub <- pt_ib[Basin == b]
  for (tt in risk_times) {
    n_risk <- sum(sub$time_years >= tt)
    risk_dt <- rbind(risk_dt, data.table(
      time = tt, Basin = b, n_risk = n_risk
    ))
  }
}
risk_dt[, Basin := factor(Basin, levels = rev(basin_order))]
risk_dt[, label := format(n_risk, big.mark = ",")]

g_risk <- ggplot(risk_dt, aes(x = time, y = Basin, label = label, color = Basin)) +
  geom_text(size = sz * 2.5 / .pt, fontface = "bold") +
  scale_color_manual(values = basin_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.5, 8), breaks = risk_times) +
  labs(x = NULL, y = NULL, title = "Number at risk") +
  theme_bw(base_size = sz) +
  theme(
    plot.title = element_text(size = sz * 2, face = "bold"),
    axis.text.x = element_text(size = sz * 3),
    axis.text.y = element_text(size = sz * 2),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 10, 10, 80)
  )

out <- file.path(base, "plot/Figures/Fig06/Fig06_diabetes_km_inter_basin.png")
png(out, width = 24, height = 18, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.75, 0.25), "npc"))))
print(g_km, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g_risk, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
cat(sprintf("Written: %s\n", out))
