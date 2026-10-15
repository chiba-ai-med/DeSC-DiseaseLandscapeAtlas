suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

out_dir <- "summary/applicability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

diseases <- c("diabetes", "dementia", "heart_failure", "sepsis")
pass_status <- c("PASS (HR=0.97)", "PASS (HR=0.99)", "FAIL (HR=1.04)", "FAIL (HR=1.06)")
names(pass_status) <- diseases

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"))

# --- 1. Energy trajectory comparison ---
traj_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_dead.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_death))
  hm_data <- traj[periods_before_death <= max_dist]

  energy_by_dist <- hm_data[, .(
    mean_energy = mean(energy, na.rm = TRUE),
    median_energy = median(energy, na.rm = TRUE),
    q25 = quantile(energy, 0.25, na.rm = TRUE),
    q75 = quantile(energy, 0.75, na.rm = TRUE)
  ), by = periods_before_death]
  setorder(energy_by_dist, periods_before_death)

  n_dead <- uniqueN(traj$kojin_id)
  energy_by_dist[, disease := d]
  energy_by_dist[, label := sprintf("%s\n%s (n=%s)",
                                     d, pass_status[d],
                                     format(n_dead, big.mark = ","))]
  energy_by_dist
})

traj_dt <- rbindlist(traj_list)
traj_dt[, label := factor(label, levels = unique(traj_dt$label))]

g_energy <- ggplot(traj_dt, aes(x = -periods_before_death)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, fill = "steelblue") +
  geom_line(aes(y = mean_energy), color = "steelblue", linewidth = 0.8) +
  geom_line(aes(y = median_energy), color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(y = mean_energy), color = "steelblue", size = 1) +
  facet_wrap(~label, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  labs(x = "Periods before death",
       y = "Energy",
       title = "A. Energy trajectory before death") +
  theme_pub +
  theme(strip.text = element_text(size = 9))

# --- 2. Basin fraction comparison ---
basin_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_dead.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_death))
  hm_data <- traj[periods_before_death <= max_dist]

  basin_by_dist <- hm_data[, .(
    frac_in_basin = mean(basin_id > 0)
  ), by = periods_before_death]
  setorder(basin_by_dist, periods_before_death)
  basin_by_dist[, disease := d]
  basin_by_dist[, label := sprintf("%s (%s)", d, pass_status[d])]
  basin_by_dist
})

basin_dt <- rbindlist(basin_list)
basin_dt[, label := factor(label, levels = unique(basin_dt$label))]

g_basin <- ggplot(basin_dt, aes(x = -periods_before_death, y = frac_in_basin)) +
  geom_col(fill = "firebrick", alpha = 0.7, width = 0.8) +
  facet_wrap(~label, nrow = 1) +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1)) +
  labs(x = "Periods before death",
       y = "In basin (%)",
       title = "B. Fraction of patients in a basin before death") +
  theme_pub +
  theme(strip.text = element_text(size = 9))

# --- 3. Heatmap comparison (top 15 states per disease) ---
hm_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_dead.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_death))
  hm_data <- traj[periods_before_death <= max_dist]

  n_at_dist <- hm_data[, .(n_patients = uniqueN(kojin_id)), by = periods_before_death]
  state_dist <- hm_data[, .(n = .N), by = .(periods_before_death, state_label)]
  state_dist <- merge(state_dist, n_at_dist, by = "periods_before_death")
  state_dist[, frac := n / n_patients]

  state_total <- state_dist[, .(total = sum(n)), by = state_label]
  setorder(state_total, -total)
  top_states <- state_total[1:min(10, nrow(state_total))]

  hm_sub <- state_dist[state_label %in% top_states$state_label]
  hm_sub[, state_label := factor(state_label,
                                   levels = rev(top_states$state_label))]
  hm_sub[, disease := d]
  hm_sub[, label := sprintf("%s (%s)", d, pass_status[d])]
  hm_sub
})

hm_dt <- rbindlist(hm_list)
hm_dt[, label := factor(label, levels = unique(hm_dt$label))]

g_heatmap <- ggplot(hm_dt, aes(x = -periods_before_death, y = state_label, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~label, nrow = 1, scales = "free_y") +
  scale_fill_viridis_c(option = "inferno", name = "Fraction",
                        begin = 0.05, end = 0.95) +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  labs(x = "Periods before death",
       y = NULL,
       title = "C. State occupancy heatmap (top 10 states)") +
  theme_pub +
  theme(strip.text = element_text(size = 9),
        axis.text.y = element_text(size = 6))

# --- Save individual panels ---
ggsave(file.path(out_dir, "panel_A_energy_trajectory.png"),
       g_energy, width = 16, height = 4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "panel_B_basin_fraction.png"),
       g_basin, width = 16, height = 3.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "panel_C_state_heatmap.png"),
       g_heatmap, width = 16, height = 6, dpi = 300, bg = "white")

cat("Written: summary/applicability/figure_representative_comparison.png\n")

# --- 4. Summary stats table ---
atlas <- fread("summary/atlas_summary.tsv")
metrics <- fread(file.path(out_dir, "disease_metrics.tsv"))

rep_tab <- metrics[disease %in% diseases,
                    .(disease, n_patients, death_rate, obs_per_patient,
                      median_dead_periods, HR, pass)]
rep_tab[, status := pass_status[disease]]
setorder(rep_tab, -pass)
fwrite(rep_tab, file.path(out_dir, "representative_diseases.tsv"), sep = "\t")
cat("Written: summary/applicability/representative_diseases.tsv\n")

cat("\nDone.\n")
