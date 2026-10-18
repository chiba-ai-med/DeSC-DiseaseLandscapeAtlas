## Applicability analysis — hospitalization version
##
## Which diseases show "low energy = high hospitalization risk" (HR < 1)?
## Logistic regression + scatter plots to explain PASS/FAIL.
##
## Usage:
##   Rscript src/applicability_analysis_hosp.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "summary/applicability_hosp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load hospitalization Cox results (whole-cohort row per disease) ---
diseases <- fread("config/disease_families.tsv")$family_id

hosp_cox <- rbindlist(lapply(diseases, function(d) {
  path <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp_specific.tsv", d)
  fread(path)
}))
hosp_cox[, pass_hosp := direction_ok == TRUE]
cat(sprintf("Hospitalization Cox: PASS %d / FAIL %d (of %d)\n",
            sum(hosp_cox$pass_hosp), sum(!hosp_cox$pass_hosp), nrow(hosp_cox)))

# --- 2. Compute per-disease metrics ---
metrics <- rbindlist(lapply(diseases, function(d) {
  pt_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", d)
  pm_path <- sprintf("output/%s/all_comorbid/cohort/person_map.tsv", d)
  traj_path <- sprintf("plot/%s/trajectories_hospitalized.tsv", d)

  pt <- fread(pt_path)
  n_patients <- nrow(pt)
  death_rate <- mean(pt$death == "dead")
  hosp_rate <- mean(pt$hospitalized == "yes")
  mean_age <- mean(pt$age)

  n_obs <- as.numeric(system(sprintf("tail -n +2 '%s' | wc -l", pm_path), intern = TRUE))
  obs_per_patient <- n_obs / n_patients

  traj <- fread(traj_path)
  median_hosp_periods <- median(traj[, .N, by = kojin_id]$N)
  n_hospitalized <- uniqueN(traj$kojin_id)

  data.table(
    disease = d,
    n_patients = n_patients,
    death_rate = death_rate,
    hosp_rate = hosp_rate,
    mean_age = mean_age,
    obs_per_patient = obs_per_patient,
    median_hosp_periods = median_hosp_periods,
    n_hospitalized = n_hospitalized
  )
}))

# Merge with hospitalization Cox
metrics <- merge(metrics,
                 hosp_cox[, .(disease_family, HR_hosp = HR, pvalue_hosp = pvalue,
                              pass_hosp, n_events, event_rate)],
                 by.x = "disease", by.y = "disease_family")

# Also merge death Cox for comparison
atlas <- fread("summary/atlas_summary.tsv")
metrics <- merge(metrics,
                 atlas[, .(disease_family, HR_death = HR, pass_death = direction_ok)],
                 by.x = "disease", by.y = "disease_family")

cat(sprintf("Death PASS: %d, FAIL: %d\n", sum(metrics$pass_death), sum(!metrics$pass_death)))
cat(sprintf("Hosp  PASS: %d, FAIL: %d\n", sum(metrics$pass_hosp), sum(!metrics$pass_hosp)))

# --- 3. Scatter plots ---
theme_pub <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

# 3a. Hosp rate vs Obs per patient
g1 <- ggplot(metrics, aes(x = obs_per_patient, y = hosp_rate,
                            color = pass_hosp, label = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Hosp Cox") +
  labs(title = "ELA applicability (hospitalization): hosp rate vs observation density",
       x = "Mean observation windows per patient",
       y = "Cohort hospitalization rate") +
  theme_pub

ggsave(file.path(out_dir, "scatter_hosprate_vs_obsperpt.png"),
       g1, width = 12, height = 8, dpi = 200, bg = "white")

# 3b. Obs per patient vs HR_hosp
g2 <- ggplot(metrics, aes(x = obs_per_patient, y = HR_hosp,
                            color = pass_hosp, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Hosp Cox") +
  labs(title = "Hospitalization HR vs observation density",
       subtitle = "HR < 1 means low energy = high hospitalization risk",
       x = "Mean observation windows per patient",
       y = "Hazard Ratio (energy, hospitalization)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_hosp_vs_obsperpt.png"),
       g2, width = 12, height = 8, dpi = 200, bg = "white")

# 3c. Death HR vs Hosp HR (comparison)
g3 <- ggplot(metrics, aes(x = HR_death, y = HR_hosp, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(aes(color = interaction(pass_death, pass_hosp)), size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 25, segment.alpha = 0.4) +
  scale_color_manual(
    values = c("TRUE.TRUE" = "steelblue", "TRUE.FALSE" = "darkorange",
               "FALSE.TRUE" = "purple", "FALSE.FALSE" = "firebrick"),
    labels = c("TRUE.TRUE" = "Both PASS", "TRUE.FALSE" = "Death PASS / Hosp FAIL",
               "FALSE.TRUE" = "Death FAIL / Hosp PASS", "FALSE.FALSE" = "Both FAIL"),
    name = NULL) +
  labs(title = "Death HR vs Hospitalization HR",
       x = "HR (energy, mortality)", y = "HR (energy, hospitalization)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_death_vs_hosp.png"),
       g3, width = 12, height = 8, dpi = 200, bg = "white")

# 3d. Mean age vs HR_hosp
g4 <- ggplot(metrics, aes(x = mean_age, y = HR_hosp,
                            color = pass_hosp, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Hosp Cox") +
  labs(title = "Hospitalization HR vs mean cohort age",
       x = "Mean age (years)", y = "Hazard Ratio (energy, hospitalization)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_hosp_vs_age.png"),
       g4, width = 12, height = 8, dpi = 200, bg = "white")

# --- 4. Logistic regression: predict hospitalization PASS/FAIL ---
metrics[, log_obs_per_pt := log(obs_per_patient)]

fit1 <- glm(pass_hosp ~ hosp_rate + log_obs_per_pt, data = metrics, family = binomial)
fit2 <- glm(pass_hosp ~ death_rate + log_obs_per_pt, data = metrics, family = binomial)
fit3 <- glm(pass_hosp ~ hosp_rate + log_obs_per_pt + mean_age, data = metrics, family = binomial)
fit4 <- glm(pass_hosp ~ obs_per_patient, data = metrics, family = binomial)

cat("\n=== Logistic: pass_hosp ~ hosp_rate + log(obs_per_patient) ===\n")
print(summary(fit1))
cat("\n=== Logistic: pass_hosp ~ death_rate + log(obs_per_patient) ===\n")
print(summary(fit2))
cat("\n=== Logistic: pass_hosp ~ hosp_rate + log(obs_per_patient) + mean_age ===\n")
print(summary(fit3))
cat("\n=== Logistic: pass_hosp ~ obs_per_patient ===\n")
print(summary(fit4))

# --- 5. Save metrics table ---
fwrite(metrics, file.path(out_dir, "disease_metrics_hosp.tsv"), sep = "\t")

# --- 6. Cross-tabulation: death vs hosp PASS/FAIL ---
cross_tab <- metrics[, .N, by = .(pass_death, pass_hosp)]
cat("\n=== Death vs Hospitalization PASS/FAIL cross-tabulation ===\n")
print(dcast(cross_tab, pass_death ~ pass_hosp, value.var = "N", fill = 0))

# --- 7. FAIL disease summary ---
fail_tab <- metrics[pass_hosp == FALSE,
  .(disease, n_patients, hosp_rate, death_rate, obs_per_patient,
    mean_age, HR_hosp, HR_death, pass_death)]
setorder(fail_tab, HR_hosp)
cat("\n=== Hospitalization FAIL diseases ===\n")
print(fail_tab)

cat("\nDone.\n")
