## Cox Energy Sanity Check per Basin — hospitalization endpoint
##
## Runs Cox (hospitalization ~ E_z + age + sex) separately for each basin.
## Uses disease-specific hospitalization (person_hospitalization_specific.tsv).
##
## Usage:
##   Rscript src/cox_per_basin_hosp.R <disease> <patient_scores_tsv> \
##     <hosp_tsv> <tekiyo_csv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
disease     <- args[1]
scores_path <- args[2]
hosp_path   <- args[3]
tekiyo_path <- args[4]
out_path    <- args[5]

cat(sprintf("=== Cox per basin (hospitalization): %s ===\n", disease))

pt <- fread(scores_path)
cat(sprintf("  %d patients\n", nrow(pt)))

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

pt[, age := start_y - parse_ym(birth_ym)]
pt[, sex_f := factor(ifelse(sex == "male", 1, 2))]

pt <- pt[is.finite(time_years) & time_years > 0 &
           is.finite(energy) & is.finite(age)]

basins <- sort(unique(pt$basin_id))
basins <- basins[basins != 0]  # exclude isolated states (basin_id=0)
n_basins <- length(basins)
cat(sprintf("  %d basins: %s\n", n_basins, paste(basins, collapse = ", ")))

results <- list()

for (b in basins) {
  sub <- pt[basin_id == b]
  n_sub <- nrow(sub)
  n_events <- sum(sub$event_hosp)
  cat(sprintf("\n  Basin %d: %d patients, %d hospitalized (%.1f%%)\n",
              b, n_sub, n_events, 100 * n_events / n_sub))

  if (n_events >= 10 && n_sub >= 30) {
    sub[, E_z := scale(energy)[, 1]]

    sd_ez <- sd(sub$E_z, na.rm = TRUE)
    if (is.na(sd_ez) || sd_ez < 1e-10) {
      cat("    Skipped: no energy variance within basin\n")
      results[[length(results) + 1]] <- data.table(
        disease_family = disease, basin_id = b, outcome = "hospitalization",
        HR = NA, HR_lower = NA, HR_upper = NA, pvalue = NA,
        direction_ok = NA, n_patients = n_sub, n_events = n_events,
        event_rate = n_events / n_sub, note = "no_energy_variance"
      )
      next
    }

    fit <- tryCatch(
      coxph(Surv(time_years, event_hosp) ~ E_z + age + sex_f, data = sub),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      s <- summary(fit)
      hr <- s$conf.int["E_z", ]
      pv <- coef(s)["E_z", 5]
      results[[length(results) + 1]] <- data.table(
        disease_family = disease, basin_id = b, outcome = "hospitalization",
        HR = hr[1], HR_lower = hr[3], HR_upper = hr[4],
        pvalue = pv, direction_ok = (hr[1] < 1.0),
        n_patients = n_sub, n_events = n_events,
        event_rate = n_events / n_sub, note = ""
      )
      cat(sprintf("    HR(E_z) = %.4f [%.4f, %.4f], p = %.2e -> %s\n",
                  hr[1], hr[3], hr[4], pv, ifelse(hr[1] < 1, "PASS", "FAIL")))
    } else {
      cat("    Cox model failed to converge\n")
      results[[length(results) + 1]] <- data.table(
        disease_family = disease, basin_id = b, outcome = "hospitalization",
        HR = NA, HR_lower = NA, HR_upper = NA, pvalue = NA,
        direction_ok = NA, n_patients = n_sub, n_events = n_events,
        event_rate = n_events / n_sub, note = "convergence_failure"
      )
    }
  } else {
    cat(sprintf("    Skipped: too few events (%d hospitalized, %d patients)\n",
                n_events, n_sub))
    results[[length(results) + 1]] <- data.table(
      disease_family = disease, basin_id = b, outcome = "hospitalization",
      HR = NA, HR_lower = NA, HR_upper = NA, pvalue = NA,
      direction_ok = NA, n_patients = n_sub, n_events = n_events,
      event_rate = n_events / n_sub, note = "too_few_events"
    )
  }
}

out <- rbindlist(results)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(out, out_path, sep = "\t")
cat(sprintf("\nWritten: %s\n", out_path))
