## Trajectory to death: track state sequence and energy for deceased patients
##
## For each dead patient, reconstruct time-ordered sequence of Allstates
## patterns across all observation windows. Outputs trajectory data and a
## heatmap showing state occupancy by distance-to-death.
##
## Usage:
##   Rscript src/trajectory_to_death.R <disease> <base_dir> <out_dir> <window_lookup_path>
##
## window_lookup_path: preprocessed TSV with columns row_id \t window
##   (created by: awk '{split($1,a,"_"); print $2 "\t" a[1]}' row_id_number_small.txt)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
disease      <- args[1]
base_dir     <- args[2]
out_dir      <- args[3]
window_path  <- args[4]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Trajectory to death: %s ===\n", disease))

# --- 1. Identify dead patients ---
pt <- fread(file.path(base_dir, "patient", "patient_scores.tsv"))
pt[, kojin_id := as.character(kojin_id)]
dead_pts <- pt[death == "dead", .(kojin_id)]
cat(sprintf("  Dead patients: %d\n", nrow(dead_pts)))

if (nrow(dead_pts) < 10) {
  cat("  Too few dead patients, skipping.\n")
  quit(save = "no")
}

# --- 2. Load person_map ---
pmap <- fread(file.path(base_dir, "cohort", "person_map.tsv"))
pmap[, kojin_id := as.character(kojin_id)]
pmap_dead <- pmap[kojin_id %chin% dead_pts$kojin_id]
cat(sprintf("  Observation rows for dead patients: %d\n", nrow(pmap_dead)))
rm(pmap); gc()

# --- 3. Map old_row_id -> time window via preprocessed lookup ---
cat("  Loading window lookup ...\n")
needed_oids <- unique(pmap_dead$old_row_id)
cat(sprintf("  Unique old_row_ids needed: %d\n", length(needed_oids)))

wlookup <- fread(window_path, header = FALSE, col.names = c("old_row_id", "window"),
                 colClasses = c("integer", "character"))
cat(sprintf("  Window lookup loaded: %d rows\n", nrow(wlookup)))

setkey(wlookup, old_row_id)
window_map <- wlookup[J(needed_oids), nomatch = NULL]
rm(wlookup); gc()
cat(sprintf("  Matched %d rows to time windows\n", nrow(window_map)))

# --- 4. Build pattern assignment + energy ---
cat("  Building pattern and energy assignment ...\n")
ls_dir <- file.path(base_dir, "landscaper", "7")
allstates <- as.matrix(fread(file.path(ls_dir, "Allstates.tsv"), header = FALSE))
subgraph <- fread(file.path(ls_dir, "SubGraph.tsv"), header = FALSE)[[1]]
energy_vec <- fread(file.path(ls_dir, "E.tsv"), header = FALSE)[[1]]
n_patterns <- nrow(allstates)

as_keys <- integer(n_patterns)
for (p in seq_len(n_patterns)) {
  as_keys[p] <- sum(2L^(which(allstates[p, ] > 0) - 1L))
}

pca_dir <- file.path(base_dir, "pca", "7")
S <- readMM(file.path(pca_dir, "Scores.mm"))
n_rows <- nrow(S)

S_bin <- as(S > 0, "TsparseMatrix")
ij_dt <- data.table(row = S_bin@i + 1L, col = S_bin@j + 1L)
ij_dt[, bit := 2L^(col - 1L)]
key_dt <- ij_dt[, .(key = sum(bit)), by = row]
rm(S, S_bin, ij_dt); gc()

row_keys <- integer(n_rows)
row_keys[key_dt$row] <- key_dt$key
rm(key_dt)

key_to_pid <- setNames(seq_len(n_patterns), as_keys)
row_pattern <- key_to_pid[match(row_keys, as_keys)]
row_pattern[is.na(row_pattern)] <- 0L
rm(row_keys)

# --- 5. Merge ---
pmap_dead[, pattern_id := row_pattern[new_row_id]]
n_subgraph <- length(subgraph)
pmap_dead[, basin_id := ifelse(pattern_id >= 1L & pattern_id <= n_subgraph,
                                subgraph[pattern_id], 0L)]
pmap_dead[is.na(basin_id), basin_id := 0L]
pmap_dead[, energy := energy_vec[pattern_id]]
pmap_dead <- merge(pmap_dead, window_map, by = "old_row_id", all.x = TRUE)
pmap_dead <- pmap_dead[!is.na(window) & pattern_id > 0]

pmap_dead[, window_start := {
  yr <- as.integer(substr(window, 1, 4))
  mo <- as.integer(substr(window, 6, 7))
  yr + (mo - 1L) / 12
}]
setorder(pmap_dead, kojin_id, window_start)

# Distance to death (0 = last observation)
pmap_dead[, n_periods := .N, by = kojin_id]
pmap_dead[, period_idx := seq_len(.N), by = kojin_id]
pmap_dead[, periods_before_death := n_periods - period_idx]

cat(sprintf("  %d trajectory rows for %d dead patients\n",
            nrow(pmap_dead), uniqueN(pmap_dead$kojin_id)))

# --- 6. State labels ---
state_labels <- character(n_patterns)
for (p in seq_len(n_patterns)) {
  signs <- ifelse(allstates[p, ] > 0, "+", "-")
  sg <- if (p <= n_subgraph) subgraph[p] else 0L
  state_labels[p] <- sprintf("S%d [%s] B%d", p, paste(signs, collapse = ""), sg)
}
pmap_dead[, state_label := state_labels[pattern_id]]
pmap_dead[, is_basin := basin_id > 0]

# --- 7. Output trajectory data ---
traj_path <- file.path(out_dir, "trajectories_dead.tsv")
fwrite(pmap_dead[, .(kojin_id, window, window_start, pattern_id,
                      basin_id, is_basin, energy, state_label,
                      periods_before_death)],
       traj_path, sep = "\t")
cat(sprintf("  Written: %s\n", traj_path))

# --- 8. Heatmap: state occupancy by distance to death ---
# Limit to last 20 periods before death
max_dist <- min(19L, max(pmap_dead$periods_before_death))
hm_data <- pmap_dead[periods_before_death <= max_dist]

# Count patients at each distance
n_at_dist <- hm_data[, .(n_patients = uniqueN(kojin_id)), by = periods_before_death]

# State frequency at each distance
state_dist <- hm_data[, .(n = .N), by = .(periods_before_death, pattern_id, state_label)]
state_dist <- merge(state_dist, n_at_dist, by = "periods_before_death")
state_dist[, frac := n / n_patients]

# Select top states by total frequency
state_total <- state_dist[, .(total = sum(n)), by = .(pattern_id, state_label)]
setorder(state_total, -total)
n_show <- min(25L, nrow(state_total))
top_states <- state_total[1:n_show]

# Mark basin states
basin_states <- unique(pmap_dead[basin_id > 0, .(pattern_id, basin_id)])
top_states <- merge(top_states, basin_states, by = "pattern_id", all.x = TRUE)
top_states[is.na(basin_id), basin_id := 0L]
top_states[, display_label := ifelse(
  basin_id > 0,
  sprintf("★ %s", state_label),
  state_label
)]

hm_plot <- merge(state_dist[pattern_id %in% top_states$pattern_id],
                 top_states[, .(pattern_id, display_label, basin_id)],
                 by = "pattern_id")

# Order: basin states first, then by total frequency
top_states[, sort_key := -total + ifelse(basin_id > 0, -1e9, 0)]
setorder(top_states, sort_key)
hm_plot[, display_label := factor(display_label, levels = top_states$display_label)]

g1 <- ggplot(hm_plot, aes(x = -periods_before_death, y = display_label, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", name = "Fraction\nof patients",
                        begin = 0.05, end = 0.95) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  labs(title = sprintf("%s: state occupancy before death (n=%d dead)",
                        disease, uniqueN(hm_data$kojin_id)),
       x = "Periods before death",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(
    size = 7,
    face = ifelse(grepl("^★", levels(hm_plot$display_label)), "bold", "plain")),
    plot.background = element_rect(fill = "white", color = NA))

hm_path <- file.path(out_dir, "trajectory_heatmap.png")
ggsave(hm_path, g1, width = 14, height = max(6, n_show * 0.35), dpi = 200, bg = "white")
cat(sprintf("  Written: %s\n", hm_path))

# --- 9. Energy trajectory plot ---
energy_by_dist <- hm_data[, .(
  mean_energy = mean(energy, na.rm = TRUE),
  median_energy = median(energy, na.rm = TRUE),
  q25 = quantile(energy, 0.25, na.rm = TRUE),
  q75 = quantile(energy, 0.75, na.rm = TRUE),
  frac_in_basin = mean(basin_id > 0)
), by = periods_before_death]
setorder(energy_by_dist, periods_before_death)

energy_by_dist_long <- melt(energy_by_dist,
  id.vars = "periods_before_death",
  measure.vars = c("mean_energy", "frac_in_basin"),
  variable.name = "metric")
energy_by_dist_long[, panel := ifelse(metric == "mean_energy", "Energy (mean)", "Fraction in basin")]

g2 <- ggplot(energy_by_dist, aes(x = -periods_before_death)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = mean_energy), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = median_energy), color = "steelblue", linetype = "dashed") +
  geom_point(aes(y = mean_energy), color = "steelblue", size = 1.5) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  labs(title = sprintf("%s: energy trajectory before death", disease),
       subtitle = "Mean (solid) and median (dashed) energy, ribbon = IQR",
       x = "Periods before death", y = "Energy") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

energy_path <- file.path(out_dir, "trajectory_energy.png")
ggsave(energy_path, g2, width = 12, height = 5, dpi = 200, bg = "white")

g2b <- ggplot(energy_by_dist, aes(x = -periods_before_death, y = frac_in_basin)) +
  geom_col(fill = "firebrick", alpha = 0.7, width = 0.8) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1)) +
  labs(title = sprintf("%s: fraction in basin before death", disease),
       x = "Periods before death", y = "In basin (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

basin_path <- file.path(out_dir, "trajectory_basin_frac.png")
ggsave(basin_path, g2b, width = 12, height = 4, dpi = 200, bg = "white")
cat(sprintf("  Written: %s\n", energy_path))
cat(sprintf("  Written: %s\n", basin_path))

# --- 10. Summary ---
cat(sprintf("\n  === Summary ===\n"))
cat(sprintf("  Dead patients tracked: %d\n", uniqueN(pmap_dead$kojin_id)))
cat(sprintf("  Median observation periods: %.0f\n",
            median(pmap_dead[, .N, by = kojin_id]$N)))
cat(sprintf("  Unique states visited: %d / %d\n",
            uniqueN(pmap_dead$pattern_id), n_patterns))

changed <- pmap_dead[, .(n_states = uniqueN(pattern_id)), by = kojin_id]
cat(sprintf("  Patients who changed state: %d / %d (%.1f%%)\n",
            sum(changed$n_states > 1), nrow(changed),
            100 * mean(changed$n_states > 1)))

at_death <- pmap_dead[periods_before_death == 0, .N, by = .(state_label, basin_id)]
setorder(at_death, -N)
cat("  Top 5 states at death:\n")
for (i in seq_len(min(5, nrow(at_death)))) {
  bflag <- ifelse(at_death$basin_id[i] > 0, " [BASIN]", "")
  cat(sprintf("    %s%s: %d (%.1f%%)\n",
              at_death$state_label[i], bflag, at_death$N[i],
              100 * at_death$N[i] / sum(at_death$N)))
}

cat("\nDone.\n")
