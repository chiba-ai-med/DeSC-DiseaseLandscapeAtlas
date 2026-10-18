## Summary plots across all 51 diseases — hospitalization endpoint
##
## Usage:
##   Rscript src/plot_entropy_summary_hosp.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "summary/applicability_hosp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

diseases <- fread("config/disease_families.tsv")$family_id

hosp_cox <- rbindlist(lapply(diseases, function(d) {
  path <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp_specific.tsv", d)
  fread(path)
}))
hosp_cox[, pass := direction_ok == TRUE]

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

# --- 1. Collect entropy data from all diseases ---
ent_all <- rbindlist(lapply(diseases, function(d) {
  f <- sprintf("plot/%s/traj_hosp_entropy.tsv", d)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, disease := d]
  dt
}))

ent_all <- merge(ent_all, hosp_cox[, .(disease_family, pass)],
                 by.x = "disease", by.y = "disease_family")
ent_all[, status := ifelse(pass, "PASS", "FAIL")]

# --- Plot 1a: All entropy trajectories overlaid ---
g1 <- ggplot(ent_all, aes(x = -periods_before_hosp, y = entropy_norm,
                           group = disease, color = status)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Hosp Cox") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 2),
    labels = function(x) ifelse(x == 0, "hosp", as.character(-x))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "State entropy trajectory before hospitalization (all 51 diseases)",
       subtitle = "Normalized Shannon entropy across 128 states",
       x = "Periods before hospitalization",
       y = "Normalized entropy") +
  theme_pub

ggsave(file.path(out_dir, "entropy_trajectories_all.png"),
       g1, width = 12, height = 6, dpi = 200, bg = "white")
cat("Written: entropy_trajectories_all.png\n")

# --- Plot 1b: Mean trajectory by PASS/FAIL ---
ent_mean <- ent_all[, .(
  mean_entropy = mean(entropy_norm, na.rm = TRUE),
  sd_entropy = sd(entropy_norm, na.rm = TRUE),
  n = .N
), by = .(periods_before_hosp, status)]

g1b <- ggplot(ent_mean, aes(x = -periods_before_hosp, color = status,
                              fill = status)) +
  geom_ribbon(aes(ymin = mean_entropy - sd_entropy,
                   ymax = mean_entropy + sd_entropy), alpha = 0.15, color = NA) +
  geom_line(aes(y = mean_entropy), linewidth = 1) +
  geom_point(aes(y = mean_entropy), size = 2) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Hosp Cox") +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Hosp Cox") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 2),
    labels = function(x) ifelse(x == 0, "hosp", as.character(-x))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Mean state entropy trajectory: PASS vs FAIL diseases (hospitalization)",
       subtitle = "Mean +/- SD across diseases",
       x = "Periods before hospitalization",
       y = "Normalized entropy") +
  theme_pub

ggsave(file.path(out_dir, "entropy_trajectory_pass_vs_fail.png"),
       g1b, width = 12, height = 6, dpi = 200, bg = "white")
cat("Written: entropy_trajectory_pass_vs_fail.png\n")

# --- 2. Compute entropy drop and non-basin fraction ---
summary_dt <- rbindlist(lapply(diseases, function(d) {
  f <- sprintf("plot/%s/traj_hosp_entropy.tsv", d)
  traj_f <- sprintf("plot/%s/trajectories_hospitalized.tsv", d)
  if (!file.exists(f) || !file.exists(traj_f)) return(NULL)

  ent <- fread(f)
  traj <- fread(traj_f)

  ent_at_hosp <- ent[periods_before_hosp == 0, entropy_norm]
  ent_early <- ent[periods_before_hosp == max(periods_before_hosp), entropy_norm]
  ent_drop <- ent_early - ent_at_hosp

  n_states_hosp <- ent[periods_before_hosp == 0, n_states]
  n_states_early <- ent[periods_before_hosp == max(periods_before_hosp), n_states]

  nonbasin_at_hosp <- mean(traj[periods_before_hosp == 0, basin_id == 0])
  nonbasin_overall <- mean(traj$basin_id == 0)

  n_hospitalized <- uniqueN(traj$kojin_id)

  data.table(
    disease = d,
    entropy_at_hosp = ent_at_hosp,
    entropy_early = ent_early,
    entropy_drop = ent_drop,
    n_states_hosp = n_states_hosp,
    n_states_early = n_states_early,
    nonbasin_at_hosp = nonbasin_at_hosp,
    nonbasin_overall = nonbasin_overall,
    n_hospitalized = n_hospitalized
  )
}))

summary_dt <- merge(summary_dt, hosp_cox[, .(disease_family, pass, HR)],
                    by.x = "disease", by.y = "disease_family")
summary_dt[, status := ifelse(pass, "PASS", "FAIL")]

metrics <- fread(file.path(out_dir, "disease_metrics_hosp.tsv"))
summary_dt <- merge(summary_dt, metrics[, .(disease, hosp_rate, obs_per_patient)],
                    by = "disease")

# --- Plot 2: Entropy drop vs hosp rate ---
g2 <- ggplot(summary_dt, aes(x = hosp_rate, y = entropy_drop,
                               color = status, label = disease)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Hosp Cox") +
  labs(title = "Entropy drop before hospitalization vs cohort hospitalization rate",
       subtitle = "Entropy drop = entropy(earliest) - entropy(hosp). Positive = concentration before hosp",
       x = "Cohort hospitalization rate",
       y = "Entropy drop (early - hosp)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_entropy_drop_vs_hosprate.png"),
       g2, width = 12, height = 8, dpi = 200, bg = "white")
cat("Written: scatter_entropy_drop_vs_hosprate.png\n")

# --- Plot 3: Non-basin fraction at hosp (bar chart) ---
setorder(summary_dt, -nonbasin_at_hosp)
summary_dt[, disease_f := factor(disease, levels = disease)]

g3 <- ggplot(summary_dt, aes(x = disease_f, y = nonbasin_at_hosp * 100,
                               fill = status)) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Hosp Cox") +
  coord_flip() +
  labs(title = "Non-basin fraction at hospitalization (all 51 diseases)",
       subtitle = "Lower = patients are mostly in basins at time of hospitalization",
       x = NULL, y = "Non-basin fraction at hospitalization (%)") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(out_dir, "barplot_nonbasin_at_hosp.png"),
       g3, width = 10, height = 12, dpi = 200, bg = "white")
cat("Written: barplot_nonbasin_at_hosp.png\n")

# --- Plot 4: Entropy at hosp vs obs_per_patient ---
g4 <- ggplot(summary_dt, aes(x = obs_per_patient, y = entropy_at_hosp,
                               color = status, label = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Hosp Cox") +
  labs(title = "State entropy at hospitalization vs observation density",
       x = "Mean observation windows per patient",
       y = "Normalized entropy at hospitalization") +
  theme_pub

ggsave(file.path(out_dir, "scatter_entropy_vs_obsperpt.png"),
       g4, width = 12, height = 8, dpi = 200, bg = "white")
cat("Written: scatter_entropy_vs_obsperpt.png\n")

# --- Plot 5: Number of states at hosp (bar chart) ---
setorder(summary_dt, -n_states_hosp)
summary_dt[, disease_f2 := factor(disease, levels = disease)]

g5 <- ggplot(summary_dt, aes(x = disease_f2, y = n_states_hosp, fill = status)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = 128, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Hosp Cox") +
  coord_flip() +
  labs(title = "Number of distinct states occupied at hospitalization (of 128 total)",
       x = NULL, y = "# distinct states") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(out_dir, "barplot_n_states_at_hosp.png"),
       g5, width = 10, height = 12, dpi = 200, bg = "white")
cat("Written: barplot_n_states_at_hosp.png\n")

# --- Save summary table ---
fwrite(summary_dt, file.path(out_dir, "entropy_summary.tsv"), sep = "\t")
cat("Written: entropy_summary.tsv\n")

cat("\nDone.\n")
