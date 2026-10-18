## Compare risk scores: Energy vs CCI vs Elixhauser
##
## Computes Harrell's C-statistic from Cox regression for each score.
##
## Usage:
##   Rscript src/compare_risk_scores.R <disease> <patient_scores> \
##       <comorbidity_tsv> <tekiyo_csv> <hosp_tsv> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
disease        <- args[1]
patient_scores <- args[2]
comorbidity_tsv <- args[3]
tekiyo_csv     <- args[4]
hosp_tsv       <- args[5]
out_dir        <- args[6]

cat(sprintf("=== Risk score comparison: %s ===\n", disease))

# --- Load data ---
ps <- fread(patient_scores)
ci <- fread(comorbidity_tsv)
tek <- fread(tekiyo_csv)

# Standardise column names
if ("kojin_id" %in% names(ps)) setnames(ps, "kojin_id", "person_id")
if ("kojin_id" %in% names(ci)) setnames(ci, "kojin_id", "person_id")

# Merge
dt <- merge(ps, ci, by = "person_id", all.x = TRUE)
dt[is.na(cci), cci := 0L]
dt[is.na(elixhauser), elixhauser := 0L]
dt[is.na(n_codes), n_codes := 0L]

# --- Death outcome ---
if ("death_date" %in% names(tek)) {
  tek_death <- tek[, .(person_id = kojin_id, death_date)]
} else if ("shibo_date" %in% names(tek)) {
  tek_death <- tek[, .(person_id = kojin_id, death_date = shibo_date)]
} else {
  tek_cols <- names(tek)
  death_col <- tek_cols[grepl("shibo|death", tek_cols, ignore.case = TRUE)][1]
  id_col <- tek_cols[grepl("kojin|person", tek_cols, ignore.case = TRUE)][1]
  tek_death <- tek[, .(person_id = get(id_col), death_date = get(death_col))]
}

dt <- merge(dt, tek_death, by = "person_id", all.x = TRUE)
dt[, dead := as.integer(!is.na(death_date) & death_date != "")]

if ("follow_up_years" %in% names(dt)) {
  dt[, time_death := follow_up_years]
} else if ("index_date" %in% names(dt) && "end_date" %in% names(dt)) {
  dt[, time_death := as.numeric(as.Date(end_date) - as.Date(index_date)) / 365.25]
} else {
  dt[, time_death := 5.0]
}
dt[time_death <= 0, time_death := 0.001]

# --- Hospitalization outcome ---
has_hosp <- file.exists(hosp_tsv) && file.size(hosp_tsv) > 0
if (has_hosp) {
  hosp <- fread(hosp_tsv)
  hosp_cols <- names(hosp)
  id_col_h <- hosp_cols[grepl("kojin|person", hosp_cols, ignore.case = TRUE)][1]
  hosp_persons <- unique(hosp[[id_col_h]])
  dt[, hospitalized := as.integer(person_id %in% hosp_persons)]
} else {
  dt[, hospitalized := 0L]
}

# --- Standardize scores for fair comparison ---
standardize <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
dt[, energy_z := standardize(energy)]
dt[, cci_z := standardize(cci)]
dt[, elix_z := standardize(elixhauser)]
dt[, ncode_z := standardize(n_codes)]
dt[, basin_f := factor(basin_id)]

# --- Compute C-statistic ---
get_cstat <- function(formula, data) {
  tryCatch({
    fit <- coxph(formula, data = data)
    conc <- concordance(fit)
    c(C = conc$concordance, se = sqrt(conc$var), HR = exp(coef(fit)[1]),
      pval = summary(fit)$coefficients[1, 5])
  }, error = function(e) c(C = NA, se = NA, HR = NA, pval = NA))
}

# --- Models ---
models <- list(
  list(name = "Energy",          formula = Surv(time_death, dead) ~ energy_z),
  list(name = "CCI",             formula = Surv(time_death, dead) ~ cci_z),
  list(name = "Elixhauser",      formula = Surv(time_death, dead) ~ elix_z),
  list(name = "N_codes",         formula = Surv(time_death, dead) ~ ncode_z),
  list(name = "Energy+Age+Sex",  formula = Surv(time_death, dead) ~ energy_z + age + sex),
  list(name = "CCI+Age+Sex",     formula = Surv(time_death, dead) ~ cci_z + age + sex),
  list(name = "Elix+Age+Sex",    formula = Surv(time_death, dead) ~ elix_z + age + sex),
  list(name = "Energy+CCI+Age+Sex",  formula = Surv(time_death, dead) ~ energy_z + cci_z + age + sex),
  list(name = "Energy+Elix+Age+Sex", formula = Surv(time_death, dead) ~ energy_z + elix_z + age + sex),
  list(name = "Basin+Age+Sex",       formula = Surv(time_death, dead) ~ basin_f + age + sex),
  list(name = "Basin+CCI+Age+Sex",   formula = Surv(time_death, dead) ~ basin_f + cci_z + age + sex),
  list(name = "Energy+Basin+CCI+Age+Sex", formula = Surv(time_death, dead) ~ energy_z + basin_f + cci_z + age + sex)
)

# Death endpoint
cat("  --- Death ---\n")
death_results <- rbindlist(lapply(models, function(m) {
  res <- get_cstat(m$formula, dt)
  cat(sprintf("    %-25s C=%.4f HR=%.4f\n", m$name, res["C"], res["HR"]))
  data.table(disease = disease, endpoint = "death", model = m$name,
             C = res["C"], se = res["se"], HR = res["HR"], pval = res["pval"])
}))

# Hospitalization endpoint
if (has_hosp && sum(dt$hospitalized) >= 20) {
  cat("  --- Hospitalization ---\n")
  models_hosp <- list(
    list(name = "Energy",          formula = Surv(time_death, hospitalized) ~ energy_z),
    list(name = "CCI",             formula = Surv(time_death, hospitalized) ~ cci_z),
    list(name = "Elixhauser",      formula = Surv(time_death, hospitalized) ~ elix_z),
    list(name = "N_codes",         formula = Surv(time_death, hospitalized) ~ ncode_z),
    list(name = "Energy+Age+Sex",  formula = Surv(time_death, hospitalized) ~ energy_z + age + sex),
    list(name = "CCI+Age+Sex",     formula = Surv(time_death, hospitalized) ~ cci_z + age + sex),
    list(name = "Elix+Age+Sex",    formula = Surv(time_death, hospitalized) ~ elix_z + age + sex),
    list(name = "Energy+CCI+Age+Sex",  formula = Surv(time_death, hospitalized) ~ energy_z + cci_z + age + sex),
    list(name = "Energy+Elix+Age+Sex", formula = Surv(time_death, hospitalized) ~ energy_z + elix_z + age + sex),
    list(name = "Basin+Age+Sex",       formula = Surv(time_death, hospitalized) ~ basin_f + age + sex),
    list(name = "Basin+CCI+Age+Sex",   formula = Surv(time_death, hospitalized) ~ basin_f + cci_z + age + sex),
    list(name = "Energy+Basin+CCI+Age+Sex", formula = Surv(time_death, hospitalized) ~ energy_z + basin_f + cci_z + age + sex)
  )
  hosp_results <- rbindlist(lapply(models_hosp, function(m) {
    res <- get_cstat(m$formula, dt)
    cat(sprintf("    %-25s C=%.4f HR=%.4f\n", m$name, res["C"], res["HR"]))
    data.table(disease = disease, endpoint = "hospitalization", model = m$name,
               C = res["C"], se = res["se"], HR = res["HR"], pval = res["pval"])
  }))
  results <- rbind(death_results, hosp_results)
} else {
  cat("  Skipping hospitalization (too few events)\n")
  results <- death_results
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "risk_score_comparison.tsv")
fwrite(results, out_path, sep = "\t")
cat(sprintf("  Written: %s\n", out_path))
cat("Done.\n")
