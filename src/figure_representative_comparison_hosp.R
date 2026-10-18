suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

out_dir <- "summary/applicability_hosp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hosp_cox <- rbindlist(lapply(
  c("diabetes", "dementia", "heart_failure", "sepsis"),
  function(d) {
    path <- sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_hosp_specific.tsv", d)
    fread(path)
  }
))

diseases <- c("diabetes", "dementia", "heart_failure", "sepsis")
pass_status <- sprintf("%s (HR=%.2f)", ifelse(hosp_cox$direction_ok, "PASS", "FAIL"), hosp_cox$HR)
names(pass_status) <- diseases

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"))

# --- 1. Energy trajectory comparison ---
traj_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_hospitalized.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_hosp))
  hm_data <- traj[periods_before_hosp <= max_dist]

  energy_by_dist <- hm_data[, .(
    mean_energy = mean(energy, na.rm = TRUE),
    median_energy = median(energy, na.rm = TRUE),
    q25 = quantile(energy, 0.25, na.rm = TRUE),
    q75 = quantile(energy, 0.75, na.rm = TRUE)
  ), by = periods_before_hosp]
  setorder(energy_by_dist, periods_before_hosp)

  n_hosp <- uniqueN(traj$kojin_id)
  energy_by_dist[, disease := d]
  energy_by_dist[, label := sprintf("%s\n%s (n=%s)",
                                     d, pass_status[d],
                                     format(n_hosp, big.mark = ","))]
  energy_by_dist
})

traj_dt <- rbindlist(traj_list)
traj_dt[, label := factor(label, levels = unique(traj_dt$label))]

g_energy <- ggplot(traj_dt, aes(x = -periods_before_hosp)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, fill = "steelblue") +
  geom_line(aes(y = mean_energy), color = "steelblue", linewidth = 0.8) +
  geom_line(aes(y = median_energy), color = "steelblue", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(y = mean_energy), color = "steelblue", size = 1) +
  facet_wrap(~label, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "hosp", as.character(-x))) +
  labs(x = "Periods before hospitalization",
       y = "Energy",
       title = "A. Energy trajectory before hospitalization") +
  theme_pub +
  theme(strip.text = element_text(size = 9))

# --- 2. Basin fraction comparison ---
basin_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_hospitalized.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_hosp))
  hm_data <- traj[periods_before_hosp <= max_dist]

  basin_by_dist <- hm_data[, .(
    frac_in_basin = mean(basin_id > 0)
  ), by = periods_before_hosp]
  setorder(basin_by_dist, periods_before_hosp)
  basin_by_dist[, disease := d]
  basin_by_dist[, label := sprintf("%s (%s)", d, pass_status[d])]
  basin_by_dist
})

basin_dt <- rbindlist(basin_list)
basin_dt[, label := factor(label, levels = unique(basin_dt$label))]

g_basin <- ggplot(basin_dt, aes(x = -periods_before_hosp, y = frac_in_basin)) +
  geom_col(fill = "firebrick", alpha = 0.7, width = 0.8) +
  facet_wrap(~label, nrow = 1) +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "hosp", as.character(-x))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1)) +
  labs(x = "Periods before hospitalization",
       y = "In basin (%)",
       title = "B. Fraction of patients in a basin before hospitalization") +
  theme_pub +
  theme(strip.text = element_text(size = 9))

# --- 3. Heatmap comparison (top 10 states per disease) ---
hm_list <- lapply(seq_along(diseases), function(i) {
  d <- diseases[i]
  traj <- fread(sprintf("plot/%s/trajectories_hospitalized.tsv", d))
  max_dist <- min(19L, max(traj$periods_before_hosp))
  hm_data <- traj[periods_before_hosp <= max_dist]

  n_at_dist <- hm_data[, .(n_patients = uniqueN(kojin_id)), by = periods_before_hosp]
  state_dist <- hm_data[, .(n = .N), by = .(periods_before_hosp, state_label)]
  state_dist <- merge(state_dist, n_at_dist, by = "periods_before_hosp")
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

g_heatmap <- ggplot(hm_dt, aes(x = -periods_before_hosp, y = state_label, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~label, nrow = 1, scales = "free_y") +
  scale_fill_viridis_c(option = "inferno", name = "Fraction",
                        begin = 0.05, end = 0.95) +
  scale_x_continuous(
    breaks = seq(-19, 0, by = 4),
    labels = function(x) ifelse(x == 0, "hosp", as.character(-x))) +
  labs(x = "Periods before hospitalization",
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

cat("Written: panel_A/B/C\n")
cat("\nDone.\n")
