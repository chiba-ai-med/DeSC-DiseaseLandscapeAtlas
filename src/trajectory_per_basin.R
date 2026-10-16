## Per-basin trajectory analysis + basin transition Sankey
##
## For each disease with 2+ basins:
##   1. Energy trajectory per basin (faceted)
##   2. Basin composition over time (stacked bar)
##   3. Basin transition Sankey (consecutive periods)
##
## Usage:
##   Rscript src/trajectory_per_basin.R <disease> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

traj_path <- sprintf("plot/%s/trajectories_dead.tsv", disease)
if (!file.exists(traj_path)) {
  cat(sprintf("No trajectory data for %s, skipping.\n", disease))
  quit(save = "no")
}

traj <- fread(traj_path)
n_basins <- length(unique(traj[basin_id > 0, basin_id]))

cat(sprintf("=== %s: %d basins, %d dead patients ===\n",
            disease, n_basins, uniqueN(traj$kojin_id)))

if (n_basins < 2) {
  cat("  Only 1 basin, skipping transition analysis.\n")
  # Still make per-basin energy trajectory (basin vs non-basin)
}

max_dist <- min(19L, max(traj$periods_before_death))
hm_data <- traj[periods_before_death <= max_dist]

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

# --- 1. Energy trajectory per basin ---
basin_label <- function(b) ifelse(b == 0, "No basin", sprintf("Basin %d", b))

energy_per_basin <- hm_data[, .(
  mean_energy = mean(energy, na.rm = TRUE),
  median_energy = median(energy, na.rm = TRUE),
  q25 = quantile(energy, 0.25, na.rm = TRUE),
  q75 = quantile(energy, 0.75, na.rm = TRUE),
  n = .N
), by = .(periods_before_death, basin_id)]

energy_per_basin[, basin_label := basin_label(basin_id)]
setorder(energy_per_basin, basin_id, periods_before_death)

# Only show basins with reasonable counts
basin_counts <- energy_per_basin[, .(total = sum(n)), by = basin_id]
keep_basins <- basin_counts[total > 100, basin_id]
epb <- energy_per_basin[basin_id %in% keep_basins]

basin_colors <- c("Basin 1" = "#E41A1C", "Basin 2" = "#377EB8",
                  "Basin 3" = "#4DAF4A", "Basin 4" = "#984EA3",
                  "No basin" = "gray50")

g1 <- ggplot(epb, aes(x = -periods_before_death, color = basin_label,
                        fill = basin_label)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.1, color = NA) +
  geom_line(aes(y = mean_energy), linewidth = 0.8) +
  geom_point(aes(y = mean_energy), size = 1.2) +
  scale_color_manual(values = basin_colors, name = NULL) +
  scale_fill_manual(values = basin_colors, name = NULL) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  labs(title = sprintf("%s: energy trajectory per basin (dead patients, n=%s)",
                        disease, format(uniqueN(hm_data$kojin_id), big.mark = ",")),
       x = "Periods before death", y = "Energy (mean, IQR ribbon)") +
  theme_pub

ggsave(file.path(out_dir, "trajectory_energy_per_basin.png"),
       g1, width = 12, height = 5, dpi = 200, bg = "white")
cat(sprintf("  Written: trajectory_energy_per_basin.png\n"))

# --- 2. Basin composition over time ---
comp <- hm_data[, .N, by = .(periods_before_death, basin_id)]
comp[, basin_label := basin_label(basin_id)]
comp_total <- comp[, .(total = sum(N)), by = periods_before_death]
comp <- merge(comp, comp_total, by = "periods_before_death")
comp[, frac := N / total]
setorder(comp, periods_before_death, basin_id)

g2 <- ggplot(comp, aes(x = -periods_before_death, y = frac,
                         fill = basin_label)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = basin_colors, name = NULL) +
  scale_x_continuous(
    breaks = seq(-max_dist, 0, by = 2),
    labels = function(x) ifelse(x == 0, "death", as.character(-x))) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(title = sprintf("%s: basin composition before death", disease),
       x = "Periods before death", y = "Fraction of patients") +
  theme_pub

ggsave(file.path(out_dir, "trajectory_basin_composition.png"),
       g2, width = 12, height = 4, dpi = 200, bg = "white")
cat(sprintf("  Written: trajectory_basin_composition.png\n"))

if (n_basins < 2) {
  cat("Done (single basin, no transition analysis).\n")
  quit(save = "no")
}

# --- 3. Basin transition analysis ---
setorder(hm_data, kojin_id, -periods_before_death)
hm_data[, basin_label := basin_label(basin_id)]

# Compute transitions between consecutive periods (earlier → later, toward death)
trans <- hm_data[, {
  if (.N >= 2) {
    data.table(
      from = basin_label[-.N],
      to = basin_label[-1],
      pbd_from = periods_before_death[-.N],
      pbd_to = periods_before_death[-1]
    )
  }
}, by = kojin_id]

# Only keep transitions between consecutive periods (diff == 1)
trans <- trans[pbd_from - pbd_to == 1L]

# Aggregate: overall transition matrix
trans_mat <- trans[, .N, by = .(from, to)]
trans_mat[, total_from := sum(N), by = from]
trans_mat[, prob := N / total_from]
setorder(trans_mat, from, -prob)

cat("\n  Transition matrix (overall):\n")
print(dcast(trans_mat, from ~ to, value.var = "prob", fill = 0))

# Save transition matrix
fwrite(trans_mat, file.path(out_dir, "basin_transition_matrix.tsv"), sep = "\t")

# --- 4. Sankey-style transition diagram ---
# Show transitions at selected time points (e.g., -15, -10, -5, -2 periods before death)
sankey_times <- sort(unique(c(
  seq(max_dist, 0, by = -5),
  2, 1, 0
)))
sankey_times <- sankey_times[sankey_times <= max_dist]
sankey_times <- sort(unique(sankey_times), decreasing = TRUE)

# For each pair of consecutive sankey times, compute transition flows
all_basins <- sort(unique(hm_data$basin_label))
n_b <- length(all_basins)

# Build node positions
node_dt <- data.table()
flow_dt <- data.table()

for (si in seq_along(sankey_times)) {
  t_val <- sankey_times[si]
  x_pos <- si

  # Count patients at this time point
  at_t <- hm_data[periods_before_death == t_val, .(n = .N), by = basin_label]
  at_t[, total := sum(n)]
  at_t[, frac := n / total]
  setorder(at_t, basin_label)

  # Build node y positions (stacked)
  y_bottom <- 0
  for (bi in seq_len(nrow(at_t))) {
    node_dt <- rbind(node_dt, data.table(
      x = x_pos, t = t_val,
      basin = at_t$basin_label[bi],
      y_bottom = y_bottom,
      y_top = y_bottom + at_t$frac[bi],
      n = at_t$n[bi],
      frac = at_t$frac[bi]
    ))
    y_bottom <- y_bottom + at_t$frac[bi]
  }

  # Compute flows to next time point
  if (si < length(sankey_times)) {
    t_next <- sankey_times[si + 1]
    # Get patients at both time points
    pts_t <- hm_data[periods_before_death == t_val, .(kojin_id, from = basin_label)]
    pts_next <- hm_data[periods_before_death == t_next, .(kojin_id, to = basin_label)]
    flows <- merge(pts_t, pts_next, by = "kojin_id")[, .N, by = .(from, to)]
    flows[, total := sum(N)]
    flows[, frac := N / total]

    flow_dt <- rbind(flow_dt, data.table(
      x_from = x_pos, x_to = x_pos + 1,
      t_from = t_val, t_to = t_next,
      from = flows$from, to = flows$to,
      n = flows$N, frac = flows$frac
    ))
  }
}

# Draw Sankey with geom_rect (nodes) + geom_polygon (flows)
# Compute flow y positions
for (si in seq_along(sankey_times)) {
  if (si >= length(sankey_times)) next
  x_l <- si
  x_r <- si + 1
  fl <- flow_dt[x_from == x_l]
  if (nrow(fl) == 0) next

  # Sort flows by from, then to
  setorder(fl, from, to)

  # For each 'from' basin, stack flows within the node's vertical space
  for (b_from in unique(fl$from)) {
    node_from <- node_dt[x == x_l & basin == b_from]
    sub_fl <- fl[from == b_from]
    sub_fl[, flow_frac := n / sum(n)]
    y_cur <- node_from$y_bottom
    for (j in seq_len(nrow(sub_fl))) {
      h <- node_from$frac * sub_fl$flow_frac[j]
      flow_dt[x_from == x_l & from == b_from & to == sub_fl$to[j],
              `:=`(y_from_bottom = y_cur, y_from_top = y_cur + h)]
      y_cur <- y_cur + h
    }
  }

  # For each 'to' basin, stack incoming flows
  for (b_to in unique(fl$to)) {
    node_to <- node_dt[x == x_r & basin == b_to]
    sub_fl <- fl[to == b_to]
    sub_fl[, flow_frac := n / sum(n)]
    y_cur <- node_to$y_bottom
    for (j in seq_len(nrow(sub_fl))) {
      h <- node_to$frac * sub_fl$flow_frac[j]
      flow_dt[x_from == x_l & to == b_to & from == sub_fl$from[j],
              `:=`(y_to_bottom = y_cur, y_to_top = y_cur + h)]
      y_cur <- y_cur + h
    }
  }
}

# Build flow polygons (curved bands)
make_band <- function(x1, x2, yb1, yt1, yb2, yt2, n_pts = 50) {
  t_seq <- seq(0, 1, length.out = n_pts)
  # Sigmoid curve
  s <- 1 / (1 + exp(-12 * (t_seq - 0.5)))
  x_top <- x1 + (x2 - x1) * t_seq
  y_top <- yt1 + (yt2 - yt1) * s
  x_bot <- rev(x_top)
  y_bot <- yb1 + (yb2 - yb1) * rev(s)
  data.table(x = c(x_top, x_bot), y = c(y_top, y_bot))
}

flow_polys <- list()
for (i in seq_len(nrow(flow_dt))) {
  r <- flow_dt[i]
  if (is.na(r$y_from_bottom)) next
  poly <- make_band(r$x_from + 0.15, r$x_to - 0.15,
                     r$y_from_bottom, r$y_from_top,
                     r$y_to_bottom, r$y_to_top)
  poly[, `:=`(id = i, from = r$from, to = r$to, n = r$n)]
  flow_polys[[i]] <- poly
}
flow_poly_dt <- rbindlist(flow_polys)

# Determine flow color (same basin = basin color, different = gray)
flow_poly_dt[, color := ifelse(from == to, from, "transition")]

# Time labels for x axis
time_labels <- ifelse(sankey_times == 0, "death", as.character(sankey_times))

g3 <- ggplot() +
  # Flow bands
  geom_polygon(data = flow_poly_dt,
               aes(x = x, y = y, group = id, fill = color),
               alpha = 0.4) +
  # Node bars
  geom_rect(data = node_dt,
            aes(xmin = x - 0.12, xmax = x + 0.12,
                ymin = y_bottom, ymax = y_top, fill = basin),
            color = "white", linewidth = 0.3) +
  # Labels on nodes
  geom_text(data = node_dt[frac > 0.05],
            aes(x = x, y = (y_bottom + y_top) / 2,
                label = sprintf("%.0f%%", frac * 100)),
            size = 2.5, fontface = "bold") +
  scale_fill_manual(
    values = c(basin_colors, "transition" = "gray70"),
    name = NULL,
    breaks = all_basins) +
  scale_x_continuous(
    breaks = seq_along(sankey_times),
    labels = time_labels) +
  scale_y_continuous(labels = NULL) +
  labs(title = sprintf("%s: basin transitions before death (n=%s)",
                        disease, format(uniqueN(hm_data$kojin_id), big.mark = ",")),
       x = "Periods before death", y = NULL) +
  theme_pub +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "trajectory_basin_sankey.png"),
       g3, width = 14, height = 5, dpi = 200, bg = "white")
cat(sprintf("  Written: trajectory_basin_sankey.png\n"))

cat("\nDone.\n")
