suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

out_dir <- "summary/applicability_hosp"

diseases <- c("diabetes", "dementia", "heart_failure", "sepsis")

hosp_cox <- rbindlist(lapply(diseases, function(d) {
  path <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp_specific.tsv", d)
  fread(path)
}))
pass_status <- sprintf("%s (HR=%.2f)", ifelse(hosp_cox$direction_ok, "PASS", "FAIL"), hosp_cox$HR)
names(pass_status) <- diseases

km_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  pt <- fread(sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", d))

  hosp <- fread(sprintf("output/%s/all_comorbid/cohort/person_hospitalization_specific.tsv", d))
  hosp[, kojin_id := as.character(kojin_id)]
  pt[, kojin_id := as.character(kojin_id)]
  pt <- merge(pt, hosp[, .(kojin_id, first_admission_ymd)],
              by = "kojin_id", all.x = TRUE)

  tek <- fread("data/tekiyo.csv",
               select = c("kojin_id", "observable_start_ym", "observable_end_ym"))
  tek[, kojin_id := as.character(kojin_id)]
  pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)

  parse_ym <- function(s) as.numeric(substr(s,1,4)) + (as.numeric(substr(s,6,7))-1)/12
  parse_ymd <- function(s) {
    y <- as.numeric(substr(s,1,4)); m <- as.numeric(substr(s,6,7))
    d2 <- ifelse(nchar(s)>=10, as.numeric(substr(s,9,10)), 15)
    y + (m-1+(d2-1)/30)/12
  }

  pt[, start_y := parse_ym(observable_start_ym)]
  pt[, end_y := parse_ym(observable_end_ym)]
  pt[, adm_y := parse_ymd(first_admission_ymd)]
  pt[, has_adm := !is.na(adm_y) & adm_y >= start_y & adm_y <= end_y]
  pt[, time_years := ifelse(has_adm, adm_y - start_y, end_y - start_y)]
  pt[, event_hosp := as.integer(has_adm)]

  pt <- pt[is.finite(time_years) & time_years > 0 & is.finite(energy)]

  brks <- quantile(pt$energy, c(0, 1/3, 2/3, 1))
  if (length(unique(brks)) < 4) return(NULL)
  pt[, energy_tertile := cut(energy, breaks = brks,
                              labels = c("low", "mid", "high"), include.lowest = TRUE)]

  sf <- survfit(Surv(time_years, event_hosp) ~ energy_tertile, data = pt)
  sdf <- data.table(
    time = sf$time, surv = sf$surv,
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
  labs(title = "B. Kaplan-Meier hospitalization-free survival by energy tertile",
       subtitle = "PASS: low energy = worse survival (HR<1);  FAIL: reversed or no separation",
       x = "Years from index date",
       y = "Hospitalization-free probability") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 9))

ggsave(file.path(out_dir, "panel_B_km_energy_tertile.png"),
       g, width = 16, height = 4.5, dpi = 300, bg = "white")

cat("Written: summary/applicability_hosp/panel_B_km_energy_tertile.png\n")
