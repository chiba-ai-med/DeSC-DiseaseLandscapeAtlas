suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

out_dir <- "summary/applicability"

diseases <- c("diabetes", "dementia", "heart_failure", "sepsis")
pass_status <- c("PASS (HR=0.97)", "PASS (HR=0.99)", "FAIL (HR=1.04)", "FAIL (HR=1.06)")
names(pass_status) <- diseases

km_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  pt <- fread(sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", d))
  pt[, dead := as.integer(death == "dead")]
  pt[, energy_tertile := cut(energy, breaks = quantile(energy, c(0, 1/3, 2/3, 1)),
                              labels = c("low", "mid", "high"), include.lowest = TRUE)]

  sf <- survfit(Surv(obs_years, dead) ~ energy_tertile, data = pt)
  sdf <- data.table(
    time = sf$time,
    surv = sf$surv,
    strata = rep(names(sf$strata), sf$strata)
  )
  sdf[, energy_tertile := sub("energy_tertile=", "", strata)]
  sdf[, disease := d]
  sdf[, label := sprintf("%s\n%s (n=%s)", d, pass_status[d],
                          format(nrow(pt), big.mark = ","))]
  sdf
})

km_dt <- rbindlist(km_list)
km_dt[, label := factor(label, levels = unique(km_dt$label))]
km_dt[, energy_tertile := factor(energy_tertile, levels = c("low", "mid", "high"))]

g <- ggplot(km_dt, aes(x = time, y = surv, color = energy_tertile)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~label, nrow = 1) +
  scale_color_manual(values = c("low" = "steelblue", "mid" = "gray50", "high" = "firebrick"),
                     name = "Energy\ntertile") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "B. Kaplan-Meier survival by energy tertile (whole cohort)",
       subtitle = "PASS: low energy = worse survival (HR<1);  FAIL: reversed or no separation",
       x = "Years from index date",
       y = "Survival probability") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 9))

ggsave(file.path(out_dir, "panel_B_km_energy_tertile.png"),
       g, width = 16, height = 4.5, dpi = 300, bg = "white")

cat("Written: summary/applicability/panel_B_km_energy_tertile.png\n")
