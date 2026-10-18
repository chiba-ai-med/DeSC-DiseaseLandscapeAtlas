## Per-disease Basin structure report
## Generates a comprehensive TSV with one row per disease
## summarizing Basin structure, outcomes, and notable codes
##
## Usage: Rscript src/disease_basin_report.R

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Per-disease Basin Report ===\n\n")

diseases <- sort(list.dirs("output", recursive = FALSE, full.names = FALSE))
diseases <- diseases[diseases != "hf_cancer"]

out_dir <- "summary/basin_exploration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reports <- list()

for (dis in diseases) {
  cat(sprintf("--- %s ---\n", dis))

  # Read patient scores
  ps_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", dis)
  if (!file.exists(ps_path)) { cat("  SKIP: no patient_scores\n"); next }
  ps <- fread(ps_path)

  # Read enrichment
  enrich_path <- sprintf("plot/%s/basin_code_enrichment.tsv", dis)
  has_enrich <- file.exists(enrich_path)
  enrich <- if (has_enrich) tryCatch(fread(enrich_path, quote = "\""), error = function(e) NULL) else NULL

  # Read comorbidity indices
  ci_path <- sprintf("output/%s/all_comorbid/comparison/comorbidity_indices.tsv", dis)
  ci <- if (file.exists(ci_path)) fread(ci_path) else NULL
  if (!is.null(ci)) {
    if ("kojin_id" %in% names(ci)) setnames(ci, "kojin_id", "person_id")
    if ("kojin_id" %in% names(ps)) setnames(ps, "kojin_id", "person_id")
    ps <- merge(ps, ci, by = "person_id", all.x = TRUE)
  }

  # Recode empty/NA basin_id to 0 (BasinNA)
  ps[is.na(basin_id) | basin_id == "", basin_id := 0L]
  ps[, basin_id := as.integer(basin_id)]

  n_total <- nrow(ps)
  n_basins <- uniqueN(ps$basin_id)
  basin_ids <- sort(unique(ps$basin_id))

  cat(sprintf("  N=%d, Basins=%d (%s)\n", n_total, n_basins, paste(basin_ids, collapse=",")))

  # Per-basin summary
  for (bid in basin_ids) {
    sub <- ps[basin_id == bid]
    n <- nrow(sub)
    pct <- 100 * n / n_total
    dr <- 100 * mean(sub$death == "dead", na.rm = TRUE)
    hr <- 100 * mean(sub$hospitalized == "yes", na.rm = TRUE)
    age_med <- median(sub$age, na.rm = TRUE)
    e_mean <- mean(sub$energy, na.rm = TRUE)
    male_pct <- 100 * mean(sub$sex == "male", na.rm = TRUE)

    cci_med <- if ("cci" %in% names(sub)) median(sub$cci, na.rm = TRUE) else NA
    elix_med <- if ("elixhauser" %in% names(sub)) median(sub$elixhauser, na.rm = TRUE) else NA
    ncode_med <- if ("n_codes" %in% names(sub)) median(sub$n_codes, na.rm = TRUE) else NA

    # Top enriched codes for this basin
    top_codes <- ""
    if (!is.null(enrich)) {
      bc <- enrich[basin_id == bid & enrichment > 1.5 & freq_in > 0.01]
      if (nrow(bc) > 0) {
        setorder(bc, -enrichment)
        top_codes <- paste(head(bc$code, 10), collapse = ",")
      }
    }

    bname <- if (bid == 0) "BasinNA" else sprintf("Basin%d", bid)
    cat(sprintf("  %s: n=%d (%.1f%%), death=%.1f%%, hosp=%.1f%%, age=%.1f, male=%.0f%%, energy=%.3f, CCI=%s, ncode=%s\n",
        bname, n, pct, dr, hr, age_med, male_pct, e_mean,
        ifelse(is.na(cci_med), "NA", sprintf("%.0f", cci_med)),
        ifelse(is.na(ncode_med), "NA", sprintf("%.0f", ncode_med))))
    if (top_codes != "") cat(sprintf("    top codes: %s\n", top_codes))

    reports[[length(reports)+1]] <- data.table(
      disease = dis, basin = bname, basin_id = bid,
      n_patients = n, pct_of_cohort = round(pct, 1),
      death_pct = round(dr, 1), hosp_pct = round(hr, 1),
      median_age = round(age_med, 1), male_pct = round(male_pct, 1),
      mean_energy = round(e_mean, 3),
      median_cci = cci_med, median_ncode = ncode_med, median_elix = elix_med,
      top_enriched_codes = top_codes
    )
  }

  # Inter-basin comparison
  if (n_basins >= 2) {
    bstats <- ps[, .(dr = mean(death == "dead", na.rm = TRUE)), by = basin_id]
    bstats <- bstats[!is.na(dr)]
    if (nrow(bstats) < 2) { cat("  Only 1 basin with data\n"); next }
    gap <- max(bstats$dr) - min(bstats$dr)
    worst <- bstats[which.max(dr)]$basin_id
    best <- bstats[which.min(dr)]$basin_id
    wname <- if (is.na(worst) || worst == 0) "BasinNA" else sprintf("Basin%d", worst)
    bname <- if (is.na(best) || best == 0) "BasinNA" else sprintf("Basin%d", best)
    cat(sprintf("  Mortality gap: %.1f%% (%s=%.1f%% vs %s=%.1f%%)\n",
        gap*100, wname, max(bstats$dr)*100, bname, min(bstats$dr)*100))
  }
  cat("\n")
}

report_dt <- rbindlist(reports)
fwrite(report_dt, file.path(out_dir, "per_disease_basin_report.tsv"), sep = "\t")
cat(sprintf("\nWritten: %s/per_disease_basin_report.tsv (%d rows)\n", out_dir, nrow(report_dt)))
cat("Done.\n")
