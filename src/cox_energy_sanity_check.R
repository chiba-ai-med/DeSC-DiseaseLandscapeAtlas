## Cox Energy Sanity Check
##
## Fits Cox proportional hazards models:
##   mortality       ~ E_z + age + sex
##   hospitalization ~ E_z + age + sex
##
## where E_z is the within-cohort Z-standardized energy.
##
## Reports HR per +1 SD Energy and flags whether the expected direction
## (HR < 1, low energy = high risk) is observed.
##
## Usage:
##   Rscript src/cox_energy_sanity_check.R \
##     <person_landscape_tsv> <hospitalization_tsv> \
##     <output_tsv> <output_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)

person_path <- args[1]
hosp_path   <- args[2]
out_tsv     <- args[3]
out_dir     <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Cox Energy Sanity Check ===\n")

# --- Load person-landscape table ---
pt <- fread(person_path)
cat(sprintf("  n = %d\n", nrow(pt)))

# --- Z-standardize energy within cohort ---
pt[, E_z := scale(energy)[, 1]]

# --- Prepare covariates ---
pt[, sex_f := factor(sex)]

# --- Filter valid observations ---
pt_valid <- pt[is.finite(time_years) & time_years > 0 & is.finite(E_z)]
cat(sprintf("  valid for Cox: n = %d\n", nrow(pt_valid)))

results <- list()

# -------------------------------------------------------------------
# Model M: mortality ~ E_z + age + sex
# -------------------------------------------------------------------
cat("\n--- Mortality ---\n")
if (sum(pt_valid$event_death, na.rm = TRUE) >= 10) {
  fit_m <- coxph(Surv(time_years, event_death) ~ E_z + age_at_baseline + sex_f,
                 data = pt_valid)
  s_m <- summary(fit_m)
  hr_m <- s_m$conf.int["E_z", ]
  pv_m <- coef(s_m)["E_z", 5]

  results$mortality <- data.table(
    outcome = "mortality",
    HR = hr_m[1],
    HR_lower = hr_m[3],
    HR_upper = hr_m[4],
    pvalue = pv_m,
    direction_ok = hr_m[1] < 1.0
  )
  cat(sprintf("  HR(E_z) = %.3f [%.3f, %.3f], p = %.2e\n",
              hr_m[1], hr_m[3], hr_m[4], pv_m))
  cat(sprintf("  direction_ok = %s\n",
              ifelse(hr_m[1] < 1.0, "TRUE", "FALSE")))

  # KM by energy tertile
  pt_valid[, E_tertile := cut(energy,
    breaks = quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("low", "mid", "high"),
    include.lowest = TRUE)]

  km <- survfit(Surv(time_years, event_death) ~ E_tertile, data = pt_valid)
  km_df <- data.frame(
    time = km$time, surv = km$surv, lower = km$lower, upper = km$upper,
    E_tertile = rep(names(km$strata), km$strata)
  )
  km_df$E_tertile <- factor(
    sub("E_tertile=", "", km_df$E_tertile),
    levels = c("low", "mid", "high"))

  g <- ggplot(km_df, aes(x = time, y = surv,
                          color = E_tertile, fill = E_tertile)) +
    geom_step(linewidth = 1.0) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    scale_color_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                   high = "#d94a3a")) +
    scale_fill_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                  high = "#d94a3a")) +
    labs(title = "Mortality KM by energy tertile (sanity check)",
         x = "Years", y = "Survival probability",
         color = "Energy", fill = "Energy") +
    theme_minimal(base_size = 14)
  ggsave(file.path(out_dir, "km_energy_mortality.png"), g,
         width = 8, height = 5, dpi = 200)
} else {
  cat("  Too few events for mortality Cox\n")
  results$mortality <- data.table(
    outcome = "mortality", HR = NA, HR_lower = NA, HR_upper = NA,
    pvalue = NA, direction_ok = NA
  )
}

# -------------------------------------------------------------------
# Model H: hospitalization ~ E_z + age + sex
# -------------------------------------------------------------------
cat("\n--- Hospitalization ---\n")
if (file.exists(hosp_path)) {
  hosp <- fread(hosp_path)
  pt_h <- merge(pt_valid, hosp, by = "kojin_id", all.x = TRUE)
  pt_h[is.na(event_hosp), event_hosp := 0L]

  if ("time_hosp_years" %in% names(pt_h)) {
    pt_h_valid <- pt_h[is.finite(time_hosp_years) & time_hosp_years > 0]
  } else {
    pt_h_valid <- pt_h[is.finite(time_years) & time_years > 0]
    pt_h_valid[, time_hosp_years := time_years]
  }

  if (sum(pt_h_valid$event_hosp, na.rm = TRUE) >= 10) {
    fit_h <- coxph(Surv(time_hosp_years, event_hosp) ~
                     E_z + age_at_baseline + sex_f,
                   data = pt_h_valid)
    s_h <- summary(fit_h)
    hr_h <- s_h$conf.int["E_z", ]
    pv_h <- coef(s_h)["E_z", 5]

    results$hospitalization <- data.table(
      outcome = "hospitalization",
      HR = hr_h[1],
      HR_lower = hr_h[3],
      HR_upper = hr_h[4],
      pvalue = pv_h,
      direction_ok = hr_h[1] < 1.0
    )
    cat(sprintf("  HR(E_z) = %.3f [%.3f, %.3f], p = %.2e\n",
                hr_h[1], hr_h[3], hr_h[4], pv_h))
    cat(sprintf("  direction_ok = %s\n",
                ifelse(hr_h[1] < 1.0, "TRUE", "FALSE")))

    # KM by energy tertile
    pt_h_valid[, E_tertile := cut(energy,
      breaks = quantile(energy, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
      labels = c("low", "mid", "high"),
      include.lowest = TRUE)]

    km_h <- survfit(Surv(time_hosp_years, event_hosp) ~ E_tertile,
                    data = pt_h_valid)
    km_h_df <- data.frame(
      time = km_h$time, surv = km_h$surv,
      lower = km_h$lower, upper = km_h$upper,
      E_tertile = rep(names(km_h$strata), km_h$strata)
    )
    km_h_df$E_tertile <- factor(
      sub("E_tertile=", "", km_h_df$E_tertile),
      levels = c("low", "mid", "high"))

    g_h <- ggplot(km_h_df, aes(x = time, y = 1 - surv,
                                color = E_tertile, fill = E_tertile)) +
      geom_step(linewidth = 1.0) +
      geom_ribbon(aes(ymin = 1 - upper, ymax = 1 - lower),
                  alpha = 0.15, color = NA) +
      scale_color_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                     high = "#d94a3a")) +
      scale_fill_manual(values = c(low = "#3a6dd9", mid = "gray50",
                                    high = "#d94a3a")) +
      labs(title = "Hospitalization cumulative incidence by energy tertile",
           x = "Years", y = "Cumulative incidence",
           color = "Energy", fill = "Energy") +
      theme_minimal(base_size = 14)
    ggsave(file.path(out_dir, "km_energy_hospitalization.png"), g_h,
           width = 8, height = 5, dpi = 200)
  } else {
    cat("  Too few hospitalization events\n")
    results$hospitalization <- data.table(
      outcome = "hospitalization", HR = NA, HR_lower = NA, HR_upper = NA,
      pvalue = NA, direction_ok = NA
    )
  }
} else {
  cat("  Hospitalization data not found, skipping\n")
  results$hospitalization <- data.table(
    outcome = "hospitalization", HR = NA, HR_lower = NA, HR_upper = NA,
    pvalue = NA, direction_ok = NA
  )
}

# --- Write results ---
result_dt <- rbindlist(results)
fwrite(result_dt, out_tsv, sep = "\t")
cat(sprintf("\nWritten: %s\n", out_tsv))

# --- Summary ---
cat("\n=== Sanity Check Summary ===\n")
for (i in seq_len(nrow(result_dt))) {
  r <- result_dt[i]
  status <- ifelse(is.na(r$direction_ok), "N/A",
                   ifelse(r$direction_ok, "PASS", "FAIL"))
  cat(sprintf("  %s: HR=%.3f [%.3f, %.3f] -> %s\n",
              r$outcome,
              ifelse(is.na(r$HR), NA, r$HR),
              ifelse(is.na(r$HR_lower), NA, r$HR_lower),
              ifelse(is.na(r$HR_upper), NA, r$HR_upper),
              status))
}
