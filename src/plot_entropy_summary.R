## Summary plots across all 51 diseases:
##   1. Entropy trajectories overlaid (PASS/FAIL colored)
##   2. Entropy drop (period 19 → death) vs death rate scatter
##   3. Non-basin fraction at death (bar chart)
##   4. Entropy at death vs obs_per_patient scatter
##
## Usage:
##   Rscript src/plot_entropy_summary.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

out_dir <- "summary/applicability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

atlas <- fread("summary/atlas_summary.tsv")
atlas[, pass := direction_ok == TRUE]

diseases <- fread("config/disease_families.tsv")$family_id

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

# --- 1. Collect entropy data from all diseases ---
ent_all <- rbindlist(lapply(diseases, function(d) {
  f <- sprintf("plot/%s/trajectory_entropy.tsv", d)
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  dt[, disease := d]
  dt
}))

ent_all <- merge(ent_all, atlas[, .(disease_family, pass)],
                 by.x = "disease", by.y = "disease_family")
ent_all[, status := ifelse(pass, "PASS", "FAIL")]

# --- Plot 1a: All entropy trajectories overlaid ---
g1 <- ggplot(ent_all, aes(x = -periods_before_death, y = entropy_norm,
                           group = disease, color = status)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Energy Cox") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "State entropy trajectory before death (all 51 diseases)",
       subtitle = "Normalized Shannon entropy across 128 states",
       x = "Periods before death",
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
), by = .(periods_before_death, status)]

g1b <- ggplot(ent_mean, aes(x = -periods_before_death, color = status,
                              fill = status)) +
  geom_ribbon(aes(ymin = mean_entropy - sd_entropy,
                   ymax = mean_entropy + sd_entropy), alpha = 0.15, color = NA) +
  geom_line(aes(y = mean_entropy), linewidth = 1) +
  geom_point(aes(y = mean_entropy), size = 2) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Energy Cox") +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Energy Cox") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Mean state entropy trajectory: PASS vs FAIL diseases",
       subtitle = "Mean +/- SD across diseases",
       x = "Periods before death",
       y = "Normalized entropy") +
  theme_pub

ggsave(file.path(out_dir, "entropy_trajectory_pass_vs_fail.png"),
       g1b, width = 12, height = 6, dpi = 200, bg = "white")
cat("Written: entropy_trajectory_pass_vs_fail.png\n")

# --- 2. Compute entropy drop and non-basin fraction ---
summary_dt <- rbindlist(lapply(diseases, function(d) {
  f <- sprintf("plot/%s/trajectory_entropy.tsv", d)
  traj_f <- sprintf("plot/%s/trajectories_dead.tsv", d)
  if (!file.exists(f) || !file.exists(traj_f)) return(NULL)

  ent <- fread(f)
  traj <- fread(traj_f)

  ent_at_death <- ent[periods_before_death == 0, entropy_norm]
  ent_early <- ent[periods_before_death == max(periods_before_death), entropy_norm]
  ent_drop <- ent_early - ent_at_death

  n_states_death <- ent[periods_before_death == 0, n_states]
  n_states_early <- ent[periods_before_death == max(periods_before_death), n_states]

  nonbasin_at_death <- mean(traj[periods_before_death == 0, basin_id == 0])
  nonbasin_overall <- mean(traj$basin_id == 0)

  n_dead <- uniqueN(traj$kojin_id)

  data.table(
    disease = d,
    entropy_at_death = ent_at_death,
    entropy_early = ent_early,
    entropy_drop = ent_drop,
    n_states_death = n_states_death,
    n_states_early = n_states_early,
    nonbasin_at_death = nonbasin_at_death,
    nonbasin_overall = nonbasin_overall,
    n_dead = n_dead
  )
}))

summary_dt <- merge(summary_dt, atlas[, .(disease_family, pass, HR)],
                    by.x = "disease", by.y = "disease_family")
summary_dt[, status := ifelse(pass, "PASS", "FAIL")]

# Load disease_metrics for death_rate and obs_per_patient
metrics <- fread(file.path(out_dir, "disease_metrics.tsv"))
summary_dt <- merge(summary_dt, metrics[, .(disease, death_rate, obs_per_patient)],
                    by = "disease")

# --- Plot 2: Entropy drop vs death rate ---
g2 <- ggplot(summary_dt, aes(x = death_rate, y = entropy_drop,
                               color = status, label = disease)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Energy Cox") +
  labs(title = "Entropy drop before death vs cohort death rate",
       subtitle = "Entropy drop = entropy(earliest) - entropy(death). Positive = concentration before death",
       x = "Cohort death rate",
       y = "Entropy drop (early - death)") +
  theme_pub

ggsave(file.path(out_dir, "scatter_entropy_drop_vs_deathrate.png"),
       g2, width = 12, height = 8, dpi = 200, bg = "white")
cat("Written: scatter_entropy_drop_vs_deathrate.png\n")

# --- Plot 3: Non-basin fraction at death (bar chart) ---
setorder(summary_dt, -nonbasin_at_death)
summary_dt[, disease_f := factor(disease, levels = disease)]

g3 <- ggplot(summary_dt, aes(x = disease_f, y = nonbasin_at_death * 100,
                               fill = status)) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Energy Cox") +
  coord_flip() +
  labs(title = "Non-basin fraction at death (all 51 diseases)",
       subtitle = "Lower = patients are mostly in basins at time of death",
       x = NULL, y = "Non-basin fraction at death (%)") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(out_dir, "barplot_nonbasin_at_death.png"),
       g3, width = 10, height = 12, dpi = 200, bg = "white")
cat("Written: barplot_nonbasin_at_death.png\n")

# --- Plot 4: Entropy at death vs obs_per_patient ---
g4 <- ggplot(summary_dt, aes(x = obs_per_patient, y = entropy_at_death,
                               color = status, label = disease)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 2.5, max.overlaps = 20, segment.alpha = 0.4) +
  scale_color_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                     name = "Energy Cox") +
  labs(title = "State entropy at death vs observation density",
       x = "Mean observation windows per patient",
       y = "Normalized entropy at death") +
  theme_pub

ggsave(file.path(out_dir, "scatter_entropy_vs_obsperpt.png"),
       g4, width = 12, height = 8, dpi = 200, bg = "white")
cat("Written: scatter_entropy_vs_obsperpt.png\n")

# --- Plot 5: Number of states at death (bar chart) ---
setorder(summary_dt, -n_states_death)
summary_dt[, disease_f2 := factor(disease, levels = disease)]

g5 <- ggplot(summary_dt, aes(x = disease_f2, y = n_states_death, fill = status)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = 128, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("PASS" = "steelblue", "FAIL" = "firebrick"),
                    name = "Energy Cox") +
  coord_flip() +
  labs(title = "Number of distinct states occupied at death (of 128 total)",
       x = NULL, y = "# distinct states") +
  theme_pub +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(out_dir, "barplot_n_states_at_death.png"),
       g5, width = 10, height = 12, dpi = 200, bg = "white")
cat("Written: barplot_n_states_at_death.png\n")

# --- Save summary table ---
fwrite(summary_dt, file.path(out_dir, "entropy_summary.tsv"), sep = "\t")
cat("Written: entropy_summary.tsv\n")

cat("\nDone.\n")
