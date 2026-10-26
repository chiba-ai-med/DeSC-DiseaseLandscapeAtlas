## Cox regression: Basin membership (BasinNA as reference) ~ mortality
## Adjusted for age + sex
## Tests whether inter-basin mortality differences survive confounding
##
## Usage: Rscript src/cox_basin_membership.R

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

cat("=== Cox: Basin membership adjusted for age + sex ===\n\n")

tekiyo_path <- "data/tekiyo.csv"
cat("Loading tekiyo ...\n")
tek <- fread(tekiyo_path,
             select = c("kojin_id", "birth_ym", "sex_code",
                         "observable_start_ym", "observable_end_ym",
                         "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
tek[, start_y := parse_ym(observable_start_ym)]
tek[, end_y := parse_ym(observable_end_ym)]
tek[, time_years := end_y - start_y]
tek[, event_death := as.integer(!is.na(shibou_flg) & shibou_flg == 1)]
tek[, age := start_y - parse_ym(birth_ym)]
tek[, sex_f := factor(sex_code)]
cat(sprintf("Before riyouka filter: %d\n", nrow(tek)))
tek <- tek[!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1]
cat(sprintf("After riyouka filter: %d\n", nrow(tek)))
tek <- tek[is.finite(time_years) & time_years > 0 & is.finite(age)]
tek <- tek[, .(kojin_id, time_years, event_death, age, sex_f)]
cat(sprintf("Tekiyo: %d patients\n\n", nrow(tek)))

diseases <- sort(list.dirs("output", recursive = FALSE, full.names = FALSE))
diseases <- diseases[diseases != "hf_cancer"]

out_dir <- "summary/basin_exploration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()

for (dis in diseases) {
  ps_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", dis)
  if (!file.exists(ps_path)) next
  ps <- fread(ps_path, select = c("kojin_id", "basin_id"))
  if (nrow(ps) == 0) next

  ps[, kojin_id := as.character(kojin_id)]
  ps[is.na(basin_id) | basin_id == "", basin_id := "0"]
  ps[, basin_id := as.character(basin_id)]

  n_basins <- uniqueN(ps$basin_id)
  if (n_basins < 2) { cat(sprintf("  %s: SKIP (1 basin)\n", dis)); next }

  dt <- merge(ps, tek, by = "kojin_id", all.x = FALSE)
  if (nrow(dt) < 100) { cat(sprintf("  %s: SKIP (too few after merge)\n", dis)); next }

  has_basinna <- "0" %in% dt$basin_id
  if (has_basinna) {
    dt[, basin_f := relevel(factor(basin_id), ref = "0")]
  } else {
    lvls <- sort(unique(dt$basin_id))
    dt[, basin_f := relevel(factor(basin_id), ref = lvls[1])]
  }

  fit <- tryCatch(
    coxph(Surv(time_years, event_death) ~ basin_f + age + sex_f, data = dt),
    error = function(e) { cat(sprintf("  %s: Cox error: %s\n", dis, e$message)); NULL }
  )
  if (is.null(fit)) next

  s <- summary(fit)
  coefs <- as.data.table(s$coefficients, keep.rownames = TRUE)
  setnames(coefs, c("term", "coef", "exp_coef", "se", "z", "p"))

  ci <- as.data.table(s$conf.int, keep.rownames = TRUE)
  setnames(ci, c("term", "HR", "neg_HR", "HR_lower", "HR_upper"))

  merged <- merge(coefs[, .(term, coef, se, z, p)],
                  ci[, .(term, HR, HR_lower, HR_upper)],
                  by = "term")

  basin_terms <- merged[grepl("^basin_f", term)]
  basin_terms[, basin_id := sub("^basin_f", "", term)]
  basin_terms[, disease := dis]
  ref_cat <- ifelse(has_basinna, "0", sort(unique(dt$basin_id))[1])
  basin_terms[, ref := ref_cat]
  basin_terms[, n_total := nrow(dt)]
  basin_terms[, n_basins := n_basins]

  cat(sprintf("  %s (ref=%s, N=%d): ", dis,
      ifelse(ref_cat == "0", "BasinNA", paste0("Basin", ref_cat)), nrow(dt)))
  for (i in seq_len(nrow(basin_terms))) {
    r <- basin_terms[i]
    sig <- ifelse(r$p < 0.001, "***", ifelse(r$p < 0.01, "**", ifelse(r$p < 0.05, "*", "")))
    cat(sprintf("Basin%s HR=%.3f(%.3f-%.3f)%s ", r$basin_id, r$HR, r$HR_lower, r$HR_upper, sig))
  }
  cat("\n")

  results[[length(results) + 1]] <- basin_terms[, .(disease, basin_id, ref, HR, HR_lower, HR_upper, p, n_total, n_basins)]
}

res_dt <- rbindlist(results, fill = TRUE)
fwrite(res_dt, file.path(out_dir, "cox_basin_membership.tsv"), sep = "\t")
cat(sprintf("\nWritten: %s/cox_basin_membership.tsv (%d rows)\n", out_dir, nrow(res_dt)))

cat("\n=== Summary: BasinNA as reference ===\n")
basinna_ref <- res_dt[ref == "0"]
if (nrow(basinna_ref) > 0) {
  cat(sprintf("  Diseases with BasinNA as ref: %d\n", uniqueN(basinna_ref$disease)))
  n_higher <- sum(basinna_ref$HR > 1 & basinna_ref$p < 0.05)
  n_lower <- sum(basinna_ref$HR < 1 & basinna_ref$p < 0.05)
  n_ns <- sum(basinna_ref$p >= 0.05)
  cat(sprintf("  Basin HR > 1 vs BasinNA (sig): %d / %d → Basins are RISKIER than BasinNA\n", n_higher, nrow(basinna_ref)))
  cat(sprintf("  Basin HR < 1 vs BasinNA (sig): %d / %d → Basins are SAFER than BasinNA\n", n_lower, nrow(basinna_ref)))
  cat(sprintf("  Not significant: %d / %d\n", n_ns, nrow(basinna_ref)))
}

cat("\nDone.\n")
