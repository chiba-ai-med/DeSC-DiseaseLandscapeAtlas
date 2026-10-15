## Cox Energy Sanity Check: Disease-specific hospitalization
##
## Uses person_hospitalization_specific.tsv (primary diagnosis matched)
## instead of all-cause hospitalization.
##
## Usage:
##   Rscript src/cox_hospitalization_specific.R <disease> <patient_scores_tsv> \
##     <hosp_specific_tsv> <tekiyo_csv> <out_tsv>

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

cat(sprintf("=== Cox disease-specific hospitalization: %s ===\n", disease))

pt <- fread(scores_path)
cat(sprintf("  %d patients\n", nrow(pt)))

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}

tek <- fread(tekiyo_path,
             select = c("kojin_id", "observable_start_ym", "observable_end_ym"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]

hosp <- fread(hosp_path)
hosp[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, hosp[, .(kojin_id, first_admission_ymd, n_admissions)],
            by = "kojin_id", all.x = TRUE, suffixes = c("", "_specific"))
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)

pt[, start_y := parse_ym(observable_start_ym)]
pt[, end_y := parse_ym(observable_end_ym)]

parse_ymd <- function(s) {
  y <- as.numeric(substr(s, 1, 4))
  m <- as.numeric(substr(s, 6, 7))
  d <- ifelse(nchar(s) >= 10, as.numeric(substr(s, 9, 10)), 15)
  y + (m - 1 + (d - 1) / 30) / 12
}
pt[, adm_y := parse_ymd(first_admission_ymd)]
pt[, has_adm := !is.na(adm_y) & adm_y >= start_y & adm_y <= end_y]
pt[, t_hosp := ifelse(has_adm, adm_y - start_y, end_y - start_y)]
pt[, event_hosp := as.integer(has_adm)]

pt_valid <- pt[is.finite(t_hosp) & t_hosp > 0 &
               is.finite(energy) & is.finite(age)]
pt_valid[, sex_f := factor(sex)]
pt_valid[, E_z := scale(energy)[, 1]]

cat(sprintf("  Valid: %d, disease-specific hosp: %d (%.1f%%)\n",
            nrow(pt_valid), sum(pt_valid$event_hosp),
            100 * mean(pt_valid$event_hosp)))

if (sum(pt_valid$event_hosp) >= 10) {
  fit <- coxph(Surv(t_hosp, event_hosp) ~ E_z + age + sex_f,
               data = pt_valid)
  s <- summary(fit)
  hr <- s$conf.int["E_z", ]
  pv <- coef(s)["E_z", 5]

  result <- data.table(
    disease_family = disease, dim = 7L, outcome = "hospitalization_specific",
    HR = hr[1], HR_lower = hr[3], HR_upper = hr[4],
    pvalue = pv, direction_ok = (hr[1] < 1.0),
    n_patients = nrow(pt_valid), n_events = sum(pt_valid$event_hosp),
    event_rate = mean(pt_valid$event_hosp)
  )
  cat(sprintf("  HR(E_z) = %.4f [%.4f, %.4f], p = %.2e -> %s\n",
              hr[1], hr[3], hr[4], pv,
              ifelse(hr[1] < 1, "PASS", "FAIL")))
} else {
  cat("  Too few events\n")
  result <- data.table(
    disease_family = disease, dim = 7L, outcome = "hospitalization_specific",
    HR = NA, HR_lower = NA, HR_upper = NA,
    pvalue = NA, direction_ok = NA,
    n_patients = nrow(pt_valid), n_events = sum(pt_valid$event_hosp),
    event_rate = mean(pt_valid$event_hosp)
  )
}

fwrite(result, out_path, sep = "\t")
cat(sprintf("Written: %s\n", out_path))
