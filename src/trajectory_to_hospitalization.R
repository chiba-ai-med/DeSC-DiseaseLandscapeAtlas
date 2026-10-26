## Trajectory to hospitalization: track state sequence and energy for hospitalized patients
##
## Same structure as trajectory_to_death.R but uses first admission as endpoint.
##
## Usage:
##   Rscript src/trajectory_to_hospitalization.R <disease> <base_dir> <out_dir> \
##     <window_lookup_path> <hosp_path>

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
hosp_path    <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== Trajectory to hospitalization: %s ===\n", disease))

# --- 1. Identify hospitalized patients ---
pt <- fread(file.path(base_dir, "patient", "patient_scores.tsv"))
pt[, kojin_id := as.character(kojin_id)]

hosp <- fread(hosp_path)
hosp[, kojin_id := as.character(kojin_id)]

hosp_pts <- pt[kojin_id %chin% hosp$kojin_id, .(kojin_id)]
hosp_pts <- merge(hosp_pts, hosp[, .(kojin_id, first_admission_ymd)], by = "kojin_id")
cat(sprintf("  Hospitalized patients: %d\n", nrow(hosp_pts)))

if (nrow(hosp_pts) < 10) {
  cat("  Too few hospitalized patients, skipping.\n")
  quit(save = "no")
}

# Compute admission window (YYYY-MM)
hosp_pts[, adm_ym := {
  y <- substr(first_admission_ymd, 1, 4)
  m <- substr(first_admission_ymd, 6, 7)
  paste0(y, "-", m)
}]

# --- 2. Load person_map ---
pmap <- fread(file.path(base_dir, "cohort", "person_map.tsv"))
pmap[, kojin_id := as.character(kojin_id)]
pmap_hosp <- pmap[kojin_id %chin% hosp_pts$kojin_id]
cat(sprintf("  Observation rows for hospitalized patients: %d\n", nrow(pmap_hosp)))
rm(pmap); gc()

# --- 3. Map old_row_id -> time window via preprocessed lookup ---
cat("  Loading window lookup ...\n")
needed_oids <- unique(pmap_hosp$old_row_id)
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
sbmap <- fread(file.path(ls_dir, "state_basin_map.tsv"))
subgraph <- sbmap$basin_id  # length 128, indexed by state_id
energy_vec <- sbmap$energy
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
pmap_hosp[, pattern_id := row_pattern[new_row_id]]
n_subgraph <- length(subgraph)
pmap_hosp[, basin_id := ifelse(pattern_id >= 1L & pattern_id <= n_subgraph,
                                subgraph[pattern_id], 0L)]
pmap_hosp[is.na(basin_id), basin_id := 0L]
pmap_hosp[, energy := energy_vec[pattern_id]]
pmap_hosp <- merge(pmap_hosp, window_map, by = "old_row_id", all.x = TRUE)
pmap_hosp <- pmap_hosp[!is.na(window) & pattern_id > 0]

# Merge admission window and filter to windows at or before admission
pmap_hosp <- merge(pmap_hosp, hosp_pts[, .(kojin_id, adm_ym)], by = "kojin_id")
pmap_hosp <- pmap_hosp[window <= adm_ym]

pmap_hosp[, window_start := {
  yr <- as.integer(substr(window, 1, 4))
  mo <- as.integer(substr(window, 6, 7))
  yr + (mo - 1L) / 12
}]
setorder(pmap_hosp, kojin_id, window_start)

# Distance to hospitalization (0 = last observation at or before admission)
pmap_hosp[, n_periods := .N, by = kojin_id]
pmap_hosp[, period_idx := seq_len(.N), by = kojin_id]
pmap_hosp[, periods_before_hosp := n_periods - period_idx]

# Drop patients with no pre-admission observations
pmap_hosp <- pmap_hosp[n_periods >= 1]

cat(sprintf("  %d trajectory rows for %d hospitalized patients\n",
            nrow(pmap_hosp), uniqueN(pmap_hosp$kojin_id)))

if (uniqueN(pmap_hosp$kojin_id) < 10) {
  cat("  Too few patients with trajectory data, skipping.\n")
  quit(save = "no")
}

# --- 6. State labels ---
state_labels <- character(n_patterns)
for (p in seq_len(n_patterns)) {
  signs <- ifelse(allstates[p, ] > 0, "+", "-")
  sg <- if (p <= n_subgraph) subgraph[p] else 0L
  state_labels[p] <- sprintf("S%d [%s] B%d", p, paste(signs, collapse = ""), sg)
}
pmap_hosp[, state_label := state_labels[pattern_id]]
pmap_hosp[, is_basin := basin_id > 0]

# --- 7. Output trajectory data ---
traj_path <- file.path(out_dir, "trajectories_hospitalized.tsv")
fwrite(pmap_hosp[, .(kojin_id, window, window_start, pattern_id,
                      basin_id, is_basin, energy, state_label,
                      periods_before_hosp)],
       traj_path, sep = "\t")
cat(sprintf("  Written: %s\n", traj_path))

# --- 8. Heatmap: state occupancy by distance to hospitalization ---
max_dist <- min(19L, max(pmap_hosp$periods_before_hosp))
hm_data <- pmap_hosp[periods_before_hosp <= max_dist]

n_at_dist <- hm_data[, .(n_patients = uniqueN(kojin_id)), by = periods_before_hosp]

state_dist <- hm_data[, .(n = .N), by = .(periods_before_hosp, pattern_id, state_label)]
state_dist <- merge(state_dist, n_at_dist, by = "periods_before_hosp")
state_dist[, frac := n / n_patients]

state_total <- state_dist[, .(total = sum(n)), by = .(pattern_id, state_label)]
setorder(state_total, -total)
n_show <- min(25L, nrow(state_total))
top_states <- state_total[1:n_show]

basin_states <- unique(pmap_hosp[basin_id > 0, .(pattern_id, basin_id)])
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

top_states[, sort_key := -total + ifelse(basin_id > 0, -1e9, 0)]
setorder(top_states, sort_key)
hm_plot[, display_label := factor(display_label, levels = top_states$display_label)]

g1 <- ggplot(hm_plot, aes(x = -periods_before_hosp, y = display_label, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", name = "Fraction\nof patients",
                        begin = 0.05, end = 0.95) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "admission", as.character(-x))) +
  labs(title = sprintf("%s: state occupancy before hospitalization (n=%d)",
                        disease, uniqueN(hm_data$kojin_id)),
       x = "Periods before first admission",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(
    size = 7,
    face = ifelse(grepl("^★", levels(hm_plot$display_label)), "bold", "plain")),
    plot.background = element_rect(fill = "white", color = NA))

hm_path <- file.path(out_dir, "traj_hosp_heatmap.png")
ggsave(hm_path, g1, width = 14, height = max(6, n_show * 0.35), dpi = 200, bg = "white")
cat(sprintf("  Written: %s\n", hm_path))

# --- 9. Energy trajectory plot ---
energy_by_dist <- hm_data[, .(
  mean_energy = mean(energy, na.rm = TRUE),
  median_energy = median(energy, na.rm = TRUE),
  q25 = quantile(energy, 0.25, na.rm = TRUE),
  q75 = quantile(energy, 0.75, na.rm = TRUE),
  frac_in_basin = mean(basin_id > 0)
), by = periods_before_hosp]
setorder(energy_by_dist, periods_before_hosp)

g2 <- ggplot(energy_by_dist, aes(x = -periods_before_hosp)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = mean_energy), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = median_energy), color = "steelblue", linetype = "dashed") +
  geom_point(aes(y = mean_energy), color = "steelblue", size = 1.5) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "admission", as.character(-x))) +
  labs(title = sprintf("%s: energy trajectory before hospitalization", disease),
       subtitle = "Mean (solid) and median (dashed) energy, ribbon = IQR",
       x = "Periods before first admission", y = "Energy") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

energy_path <- file.path(out_dir, "traj_hosp_energy.png")
ggsave(energy_path, g2, width = 12, height = 5, dpi = 200, bg = "white")

g2b <- ggplot(energy_by_dist, aes(x = -periods_before_hosp, y = frac_in_basin)) +
  geom_col(fill = "firebrick", alpha = 0.7, width = 0.8) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "admission", as.character(-x))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                     limits = c(0, 1)) +
  labs(title = sprintf("%s: fraction in basin before hospitalization", disease),
       x = "Periods before first admission", y = "In basin (%)") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA))

basin_path <- file.path(out_dir, "traj_hosp_basin_frac.png")
ggsave(basin_path, g2b, width = 12, height = 4, dpi = 200, bg = "white")
cat(sprintf("  Written: %s\n", energy_path))
cat(sprintf("  Written: %s\n", basin_path))

# --- 10. Summary ---
cat(sprintf("\n  === Summary ===\n"))
cat(sprintf("  Hospitalized patients tracked: %d\n", uniqueN(pmap_hosp$kojin_id)))
cat(sprintf("  Median observation periods: %.0f\n",
            median(pmap_hosp[, .N, by = kojin_id]$N)))
cat(sprintf("  Unique states visited: %d / %d\n",
            uniqueN(pmap_hosp$pattern_id), n_patterns))

changed <- pmap_hosp[, .(n_states = uniqueN(pattern_id)), by = kojin_id]
cat(sprintf("  Patients who changed state: %d / %d (%.1f%%)\n",
            sum(changed$n_states > 1), nrow(changed),
            100 * mean(changed$n_states > 1)))

at_hosp <- pmap_hosp[periods_before_hosp == 0, .N, by = .(state_label, basin_id)]
setorder(at_hosp, -N)
cat("  Top 5 states at hospitalization:\n")
for (i in seq_len(min(5, nrow(at_hosp)))) {
  bflag <- ifelse(at_hosp$basin_id[i] > 0, " [BASIN]", "")
  cat(sprintf("    %s%s: %d (%.1f%%)\n",
              at_hosp$state_label[i], bflag, at_hosp$N[i],
              100 * at_hosp$N[i] / sum(at_hosp$N)))
}

cat("\nDone.\n")
