## Cox model comparison: energy only vs energy + n_diagnoses
##
## Test whether energy's mortality prediction survives after controlling
## for total diagnosis count (n_diagnoses).
##
## Usage:
##   Rscript src/cox_with_ndiag.R <disease_family>

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

fid <- commandArgs(trailingOnly = TRUE)[1]

scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
ndiag_path <- sprintf("output/%s/all_comorbid/cohort/n_diagnoses.tsv", fid)
out_dir <- sprintf("summary/%s_deep_dive", fid)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pt <- fread(scores_path)
pt <- pt[is.finite(energy)]
pt[, death_bin := as.integer(death == "dead")]
pt[, sex_f := as.integer(sex == "female")]

ndiag <- fread(ndiag_path)
pt <- merge(pt, ndiag, by = "kojin_id", all.x = TRUE)
pt <- pt[!is.na(n_diagnoses)]

cat(sprintf("Patients with n_diagnoses: %d\n", nrow(pt)))
cat(sprintf("cor(energy, n_diagnoses): %.4f\n",
            cor(pt$energy, pt$n_diagnoses, use = "pairwise.complete.obs")))

# Standardize
pt[, E_z := scale(energy)]
pt[, ndiag_z := scale(n_diagnoses)]

pt <- pt[obs_years > 0]

# --- Model 1: Original (energy + age + sex) ---
m1 <- coxph(Surv(obs_years, death_bin) ~ E_z + age + sex_f, data = pt)
s1 <- summary(m1)
cat("\n=== Model 1: Energy + Age + Sex ===\n")
print(s1$coefficients)

# --- Model 2: Energy + Age + Sex + n_diagnoses ---
m2 <- coxph(Surv(obs_years, death_bin) ~ E_z + age + sex_f + ndiag_z, data = pt)
s2 <- summary(m2)
cat("\n=== Model 2: Energy + Age + Sex + n_diagnoses ===\n")
print(s2$coefficients)

# --- Model 3: n_diagnoses only (for comparison) ---
m3 <- coxph(Surv(obs_years, death_bin) ~ ndiag_z + age + sex_f, data = pt)
s3 <- summary(m3)
cat("\n=== Model 3: n_diagnoses + Age + Sex (no energy) ===\n")
print(s3$coefficients)

# --- Summary comparison ---
cat("\n=== Comparison ===\n")
cat(sprintf("Model 1 (E+age+sex):       HR(E)=%.4f, p=%.2e, concordance=%.4f\n",
            exp(s1$coefficients["E_z", "coef"]),
            s1$coefficients["E_z", "Pr(>|z|)"],
            s1$concordance["C"]))
cat(sprintf("Model 2 (E+age+sex+ndiag): HR(E)=%.4f, p=%.2e, concordance=%.4f\n",
            exp(s2$coefficients["E_z", "coef"]),
            s2$coefficients["E_z", "Pr(>|z|)"],
            s2$concordance["C"]))
cat(sprintf("Model 3 (ndiag+age+sex):   HR(ndiag)=%.4f, p=%.2e, concordance=%.4f\n",
            exp(s3$coefficients["ndiag_z", "coef"]),
            s3$coefficients["ndiag_z", "Pr(>|z|)"],
            s3$concordance["C"]))

# Energy tertile x n_diagnoses
pt[, E_tertile := cut(energy, breaks = quantile(energy, c(0, 1/3, 2/3, 1)),
                       labels = c("Low", "Mid", "High"), include.lowest = TRUE)]
cat("\n=== n_diagnoses by energy tertile ===\n")
print(pt[, .(mean_ndiag = mean(n_diagnoses), median_ndiag = median(n_diagnoses),
             sd_ndiag = sd(n_diagnoses), death_rate = mean(death_bin)),
         by = E_tertile][order(E_tertile)])

# Save results
results <- data.table(
  model = c("E+age+sex", "E+age+sex+ndiag", "ndiag+age+sex"),
  HR_energy = c(exp(s1$coefficients["E_z", "coef"]),
                exp(s2$coefficients["E_z", "coef"]),
                NA),
  p_energy = c(s1$coefficients["E_z", "Pr(>|z|)"],
               s2$coefficients["E_z", "Pr(>|z|)"],
               NA),
  HR_ndiag = c(NA,
               exp(s2$coefficients["ndiag_z", "coef"]),
               exp(s3$coefficients["ndiag_z", "coef"])),
  p_ndiag = c(NA,
              s2$coefficients["ndiag_z", "Pr(>|z|)"],
              s3$coefficients["ndiag_z", "Pr(>|z|)"]),
  concordance = c(s1$concordance["C"], s2$concordance["C"], s3$concordance["C"])
)
fwrite(results, file.path(out_dir, "cox_energy_vs_ndiag.tsv"), sep = "\t")
cat(sprintf("\nOutput: %s/cox_energy_vs_ndiag.tsv\n", out_dir))
