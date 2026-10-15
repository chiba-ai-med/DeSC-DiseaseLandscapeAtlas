suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "summary/applicability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load atlas summary (PASS/FAIL) ---
atlas <- fread("summary/atlas_summary.tsv")
atlas[, pass := direction_ok == TRUE]

# --- 2. Compute per-disease metrics ---
diseases <- fread("config/disease_families.tsv")$family_id

metrics <- rbindlist(lapply(diseases, function(d) {
  pt_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", d)
  pm_path <- sprintf("output/%s/all_comorbid/cohort/person_map.tsv", d)
  traj_path <- sprintf("plot/%s/trajectories_dead.tsv", d)

  pt <- fread(pt_path)
  n_patients <- nrow(pt)
  death_rate <- mean(pt$death == "dead")
  mean_obs_years <- mean(pt$obs_years)
  mean_age <- mean(pt$age)

  n_obs <- as.numeric(system(sprintf("tail -n +2 '%s' | wc -l", pm_path), intern = TRUE))
  obs_per_patient <- n_obs / n_patients

  traj <- fread(traj_path)
  median_dead_periods <- median(traj[, .N, by = kojin_id]$N)
  n_dead <- uniqueN(traj$kojin_id)

  data.table(
    disease = d,
    n_patients = n_patients,
    death_rate = death_rate,
    mean_obs_years = mean_obs_years,
    mean_age = mean_age,
    obs_per_patient = obs_per_patient,
    median_dead_periods = median_dead_periods,
    n_dead = n_dead
  )
}))

# Merge with atlas
metrics <- merge(metrics, atlas[, .(disease_family, HR, pvalue, pass, n_basins,
                                     largest_basin_frac)],
                 by.x = "disease", by.y = "disease_family")

cat(sprintf("PASS: %d, FAIL: %d\n", sum(metrics$pass), sum(!metrics$pass)))

# --- 3. Scatter plots ---
theme_pub <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

# 3a. Death rate vs Obs per patient
g1 <- ggplot(metrics, aes(x = obs_per_patient, y = death_rate,
                            color = pass, label = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Energy Cox\nsanity check") +
  labs(title = "ELA applicability: death rate vs observation density",
       x = "Mean observation windows per patient",
       y = "Cohort death rate") +
  theme_pub

ggsave(file.path(out_dir, "scatter_deathrate_vs_obsperpt.png"),
       g1, width = 12, height = 8, dpi = 200, bg = "white")

# 3b. Death rate vs median dead periods
g2 <- ggplot(metrics, aes(x = median_dead_periods, y = death_rate,
                            color = pass, label = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Energy Cox\nsanity check") +
  labs(title = "ELA applicability: death rate vs dead-patient observation depth",
       x = "Median observation periods (dead patients)",
       y = "Cohort death rate") +
  theme_pub

ggsave(file.path(out_dir, "scatter_deathrate_vs_deadperiods.png"),
       g2, width = 12, height = 8, dpi = 200, bg = "white")

# 3c. Obs per patient vs HR (continuous)
g3 <- ggplot(metrics, aes(x = obs_per_patient, y = HR,
                            color = pass, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Energy Cox\nsanity check") +
  labs(title = "Energy HR vs observation density",
       subtitle = "HR < 1 means low energy = high risk (expected)",
       x = "Mean observation windows per patient",
       y = "Hazard Ratio (energy, mortality)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_vs_obsperpt.png"),
       g3, width = 12, height = 8, dpi = 200, bg = "white")

# 3d. Mean age vs HR
g4 <- ggplot(metrics, aes(x = mean_age, y = HR,
                            color = pass, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Energy Cox\nsanity check") +
  labs(title = "Energy HR vs mean cohort age",
       x = "Mean age (years)", y = "Hazard Ratio (energy, mortality)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_vs_age.png"),
       g4, width = 12, height = 8, dpi = 200, bg = "white")

# 3e. n_basins vs HR
g5 <- ggplot(metrics, aes(x = n_basins, y = HR,
                            color = pass, label = disease)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(width = 0.1)) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
                     labels = c("TRUE" = "PASS", "FALSE" = "FAIL"),
                     name = "Energy Cox\nsanity check") +
  labs(title = "Energy HR vs number of basins",
       x = "Number of basins", y = "Hazard Ratio (energy, mortality)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_HR_vs_nbasins.png"),
       g5, width = 12, height = 8, dpi = 200, bg = "white")

# --- 4. Logistic regression: predict PASS/FAIL ---
metrics[, log_obs_per_pt := log(obs_per_patient)]

fit1 <- glm(pass ~ death_rate + log_obs_per_pt, data = metrics, family = binomial)
fit2 <- glm(pass ~ death_rate + median_dead_periods, data = metrics, family = binomial)
fit3 <- glm(pass ~ death_rate + log_obs_per_pt + mean_age, data = metrics, family = binomial)
fit4 <- glm(pass ~ obs_per_patient, data = metrics, family = binomial)

cat("\n=== Logistic regression: PASS ~ death_rate + log(obs_per_patient) ===\n")
print(summary(fit1))
cat("\n=== Logistic regression: PASS ~ death_rate + median_dead_periods ===\n")
print(summary(fit2))
cat("\n=== Logistic regression: PASS ~ death_rate + log(obs_per_patient) + mean_age ===\n")
print(summary(fit3))
cat("\n=== Logistic regression: PASS ~ obs_per_patient ===\n")
print(summary(fit4))

# --- 5. Save metrics table ---
fwrite(metrics, file.path(out_dir, "disease_metrics.tsv"), sep = "\t")

# --- 6. FAIL disease summary table ---
fail_tab <- metrics[pass == FALSE, .(disease, n_patients, death_rate, obs_per_patient,
                                      median_dead_periods, mean_age, HR, n_basins)]
setorder(fail_tab, HR)
cat("\n=== FAIL diseases ===\n")
print(fail_tab)

cat("\nDone.\n")
