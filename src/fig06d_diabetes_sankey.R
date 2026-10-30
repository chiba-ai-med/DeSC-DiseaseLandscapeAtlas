suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
disease <- "diabetes"

traj <- fread(file.path(base, sprintf("plot/%s/trajectories_dead.tsv", disease)))
n_basins <- length(unique(traj[basin_id > 0, basin_id]))

cat(sprintf("=== %s: %d basins, %d dead patients ===\n",
            disease, n_basins, uniqueN(traj$kojin_id)))

max_dist <- min(19L, max(traj$periods_before_death))
hm_data <- traj[periods_before_death <= max_dist]

theme_pub <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank())

basin_label_fn <- function(b) {
  ifelse(b == 0, "Non-basin", sprintf("Basin %d", b))
}

hm_data[, group := basin_label_fn(basin_id)]

basin_colors <- c("Basin 1" = "#E41A1C", "Basin 2" = "#377EB8",
                  "Basin 3" = "#4DAF4A", "Basin 4" = "#984EA3",
                  "Non-basin" = "gray50")

sankey_times <- sort(unique(c(
  seq(max_dist, 0, by = -5), 2, 1, 0
)), decreasing = TRUE)

all_groups <- sort(unique(hm_data$group))

node_dt <- data.table()
for (si in seq_along(sankey_times)) {
  t_val <- sankey_times[si]
  at_t <- hm_data[periods_before_death == t_val, .(n = .N), by = group]
  at_t[, total := sum(n)]
  at_t[, frac := n / total]
  for (g in all_groups) {
    if (!g %in% at_t$group) {
      at_t <- rbind(at_t, data.table(group = g, n = 0L, total = at_t$total[1], frac = 0))
    }
  }
  setorder(at_t, group)
  y_bottom <- 0
  for (bi in seq_len(nrow(at_t))) {
    node_dt <- rbind(node_dt, data.table(
      x = si, t = t_val,
      group = at_t$group[bi],
      y_bottom = y_bottom,
      y_top = y_bottom + at_t$frac[bi],
      n = at_t$n[bi],
      frac = at_t$frac[bi]
    ))
    y_bottom <- y_bottom + at_t$frac[bi]
  }
}

setorder(hm_data, kojin_id, -periods_before_death)
flow_dt <- data.table()
for (si in seq_len(length(sankey_times) - 1)) {
  t_from <- sankey_times[si]
  t_to <- sankey_times[si + 1]

  pts_from <- hm_data[periods_before_death == t_from, .(kojin_id, from = group)]
  pts_to <- hm_data[periods_before_death == t_to, .(kojin_id, to = group)]
  fl <- merge(pts_from, pts_to, by = "kojin_id")[, .N, by = .(from, to)]
  fl[, total := sum(N)]
  fl[, frac := N / total]

  flow_dt <- rbind(flow_dt, data.table(
    x_from = si, x_to = si + 1,
    t_from = t_from, t_to = t_to,
    from = fl$from, to = fl$to,
    n = fl$N, frac = fl$frac
  ))
}

for (si in seq_len(length(sankey_times) - 1)) {
  fl <- flow_dt[x_from == si]
  if (nrow(fl) == 0) next
  setorder(fl, from, to)

  for (b_from in unique(fl$from)) {
    node_from <- node_dt[x == si & group == b_from]
    if (nrow(node_from) == 0 || node_from$frac == 0) next
    sub_fl <- fl[from == b_from]
    sub_fl[, flow_frac := n / sum(n)]
    y_cur <- node_from$y_bottom
    for (j in seq_len(nrow(sub_fl))) {
      h <- node_from$frac * sub_fl$flow_frac[j]
      flow_dt[x_from == si & from == b_from & to == sub_fl$to[j],
              `:=`(y_from_bottom = y_cur, y_from_top = y_cur + h)]
      y_cur <- y_cur + h
    }
  }

  for (b_to in unique(fl$to)) {
    node_to <- node_dt[x == si + 1 & group == b_to]
    if (nrow(node_to) == 0 || node_to$frac == 0) next
    sub_fl <- fl[to == b_to]
    sub_fl[, flow_frac := n / sum(n)]
    y_cur <- node_to$y_bottom
    for (j in seq_len(nrow(sub_fl))) {
      h <- node_to$frac * sub_fl$flow_frac[j]
      flow_dt[x_from == si & to == b_to & from == sub_fl$from[j],
              `:=`(y_to_bottom = y_cur, y_to_top = y_cur + h)]
      y_cur <- y_cur + h
    }
  }
}

make_band <- function(x1, x2, yb1, yt1, yb2, yt2, n_pts = 50) {
  t_seq <- seq(0, 1, length.out = n_pts)
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
  if (r$y_from_top - r$y_from_bottom < 1e-6) next
  poly <- make_band(r$x_from + 0.15, r$x_to - 0.15,
                     r$y_from_bottom, r$y_from_top,
                     r$y_to_bottom, r$y_to_top)
  poly[, `:=`(id = i, from = r$from, to = r$to)]
  flow_polys[[i]] <- poly
}
flow_poly_dt <- rbindlist(flow_polys)
flow_poly_dt[, color := ifelse(from == to, from, "transition")]

time_labels <- ifelse(sankey_times == 0, "death", as.character(sankey_times))

dummy_rect <- data.table(xmin = -99, xmax = -98, ymin = -99, ymax = -98,
                         group = "Non-basin")

g2 <- ggplot() +
  geom_rect(data = dummy_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group)) +
  geom_polygon(data = flow_poly_dt,
               aes(x = x, y = y, group = id, fill = color),
               alpha = 0.4) +
  geom_rect(data = node_dt[frac > 0],
            aes(xmin = x - 0.12, xmax = x + 0.12,
                ymin = y_bottom, ymax = y_top, fill = group),
            color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(basin_colors, "transition" = "gray70"),
    name = NULL,
    breaks = c("Basin 1", "Basin 2", "Basin 3", "Basin 4", "Non-basin"),
    drop = FALSE) +
  scale_x_continuous(breaks = seq_along(sankey_times), labels = time_labels) +
  scale_y_continuous(labels = NULL) +
  coord_cartesian(xlim = c(0.5, length(sankey_times) + 0.5), ylim = c(0, 1)) +
  labs(title = sprintf("Diabetes: basin transitions before death (n=%s)",
                        format(uniqueN(hm_data$kojin_id), big.mark = ",")),
       x = "Periods before death", y = NULL) +
  theme_pub +
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 11 * 3),
        axis.title.x = element_text(size = 11 * 3))

out <- file.path(base, "plot/Figures/Fig06/Fig06_diabetes_sankey.png")
ggsave(out, g2, width = 14, height = 5 * 1.1, dpi = 200, bg = "white")
cat(sprintf("Written: %s\n", out))
