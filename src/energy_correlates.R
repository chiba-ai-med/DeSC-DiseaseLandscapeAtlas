## What does Energy correlate with?
## Compute Spearman correlations between Energy and patient metadata
## across all 51 diseases.
##
## Usage:
##   Rscript src/energy_correlates.R

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== Energy correlates analysis ===\n")

diseases <- list.dirs("output", recursive = FALSE, full.names = FALSE)
diseases <- diseases[diseases != "hf_cancer"]

results <- list()

for (dis in diseases) {
  scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", dis)
  indices_path <- sprintf("output/%s/all_comorbid/comparison/comorbidity_indices.tsv", dis)

  if (!file.exists(scores_path) || !file.exists(indices_path)) next

  dt <- fread(scores_path)
  idx <- fread(indices_path)
  dt <- merge(dt, idx, by = "kojin_id", all.x = TRUE)

  if (nrow(dt) < 100 || all(is.na(dt$energy))) next

  # Numeric variables to correlate with Energy
  dt[, sex_num := ifelse(sex == "male", 1, 0)]
  dt[, death_num := ifelse(death == "dead", 1, 0)]
  dt[, hosp_num := ifelse(hospitalized == "yes", 1, 0)]

  vars <- c("age", "sex_num", "n_codes", "cci", "elixhauser",
            "n_admissions", "obs_years", "death_num", "hosp_num",
            "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")

  cors <- sapply(vars, function(v) {
    if (!v %in% names(dt)) return(NA_real_)
    x <- dt[[v]]
    y <- dt$energy
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 50) return(NA_real_)
    cor(x[ok], y[ok], method = "spearman")
  })

  results[[dis]] <- data.table(disease = dis, variable = names(cors), rho = cors)
}

all_cors <- rbindlist(results)

# Summary: median and IQR of correlation across diseases
summary_dt <- all_cors[, .(
  median_rho = median(rho, na.rm = TRUE),
  q25 = quantile(rho, 0.25, na.rm = TRUE),
  q75 = quantile(rho, 0.75, na.rm = TRUE),
  n_diseases = sum(!is.na(rho))
), by = variable]
setorder(summary_dt, -median_rho)

cat("\n=== Spearman correlation with Energy (51 diseases, median [IQR]) ===\n")
for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  cat(sprintf("  %-15s  rho = %+.3f [%+.3f, %+.3f]  (n=%d)\n",
      r$variable, r$median_rho, r$q25, r$q75, r$n_diseases))
}

out_dir <- "summary/risk_comparison"
fwrite(all_cors, file.path(out_dir, "energy_correlates.tsv"), sep = "\t")
fwrite(summary_dt, file.path(out_dir, "energy_correlates_summary.tsv"), sep = "\t")
cat(sprintf("\nWritten to %s\n", out_dir))
cat("Done.\n")
