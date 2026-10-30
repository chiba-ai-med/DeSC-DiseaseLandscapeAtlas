## Generate case study panels for any disease
## Usage: Rscript src/fig_casestudy_panels.R <disease> <fig_num> <basin_a> <basin_b>
##   basin_a, basin_b: which basins to compare in panel a

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
disease   <- args[1]
fig_num   <- args[2]
basin_a   <- as.integer(args[3])
basin_b   <- as.integer(args[4])

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
fig_dir <- file.path(base, sprintf("plot/Figures/Fig%s", fig_num))
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

disease_label <- gsub("_", " ", disease)
disease_label <- paste0(toupper(substr(disease_label, 1, 1)),
                        substr(disease_label, 2, nchar(disease_label)))

cat(sprintf("=== %s (Fig%s) ===\n", disease_label, fig_num))

## ── Panel a: logitdiff comparison ──
tsv <- file.path(base, sprintf("plot/%s/tagcloud_scored/logitdiff/basin_code_logitdiff.tsv", disease))
dt <- fread(tsv)

ps <- fread(file.path(base, sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", disease)),
            select = "basin_id")
bpat <- ps[basin_id > 0, .(n_patients = .N), by = basin_id]

cox <- fread(file.path(base, sprintf("output/%s/all_comorbid/cox/7/energy_sanity_check_per_basin.tsv", disease)))

make_basin_label <- function(bid) {
  n <- bpat[basin_id == bid, n_patients]
  hr <- cox[basin_id == bid & outcome == "mortality", round(HR, 2)]
  sprintf("Basin %d  HR=%s (n=%s)", bid, hr, format(n, big.mark = ","))
}

top_n <- 10
truncate_label <- function(s, maxlen = 40) {
  ifelse(nchar(s) > maxlen, paste0(substr(s, 1, maxlen - 1), "…"), s)
}

ba <- dt[basin_id == basin_a & freq_in > 0.01][order(-score)][1:top_n]
bb <- dt[basin_id == basin_b & freq_in > 0.01][order(-score)][1:top_n]

ba[, label := paste0(truncate_label(en_label), " (", code, ")")]
bb[, label := paste0(truncate_label(en_label), " (", code, ")")]

lab_a <- make_basin_label(basin_a)
lab_b <- make_basin_label(basin_b)
ba[, basin := lab_a]
bb[, basin := lab_b]
ba[, label := factor(label, levels = rev(label))]
bb[, label := factor(label, levels = rev(label))]

combined <- rbindlist(list(ba, bb))
combined[, basin := factor(basin, levels = c(lab_a, lab_b))]

sz <- 16
g_a <- ggplot(combined, aes(x = label, y = score, fill = freq_in)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", freq_in * 100)),
            hjust = -0.1, size = sz / 3.5, fontface = "bold") +
  coord_flip() +
  facet_wrap(~basin, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Freq",
                       begin = 0.15, end = 0.95,
                       labels = scales::percent) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = sprintf("%s - Top 10 comorbidity codes by basin", disease_label),
       x = NULL, y = "Log-odds ratio") +
  theme_bw(base_size = sz) +
  theme(plot.title = element_text(face = "bold", size = sz),
        axis.text.y = element_text(size = sz * 0.9),
        axis.text.x = element_text(size = sz),
        axis.title.x = element_text(size = sz),
        strip.text = element_text(size = sz * 0.75, face = "bold"),
        legend.text = element_text(size = sz - 2),
        legend.title = element_text(size = sz - 2),
        plot.background = element_rect(fill = "white", color = NA))

fig_h <- (1.5 + 0.7 * top_n) * 0.8
out_a <- file.path(fig_dir, sprintf("Fig%s_a_%s_logitdiff.png", fig_num, disease))
ggsave(out_a, g_a, width = 16 * 0.8, height = fig_h, dpi = 300, bg = "white")
cat(sprintf("  Written: %s\n", out_a))

## ── Panel b: KM inter-basin ──
tekiyo_path <- file.path(base, "data/tekiyo.csv")
pt <- fread(file.path(base, sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", disease)))
tek <- fread(tekiyo_path,
             select = c("kojin_id", "observable_start_ym", "observable_end_ym",
                        "shibou_flg", "shibou_flg_riyouka"))
tek[, kojin_id := as.character(kojin_id)]
pt[, kojin_id := as.character(kojin_id)]
pt <- merge(pt, tek, by = "kojin_id", all.x = TRUE)
rm(tek)

parse_ym <- function(s) {
  as.numeric(substr(s, 1, 4)) + (as.numeric(substr(s, 6, 7)) - 1) / 12
}
pt[, time_years := parse_ym(observable_end_ym) - parse_ym(observable_start_ym)]
pt[, event_death := as.integer(!is.na(shibou_flg) & shibou_flg == 1)]
pt <- pt[!is.na(shibou_flg_riyouka) & shibou_flg_riyouka == 1]
pt <- pt[is.finite(time_years) & time_years > 0 & is.finite(energy)]

pt_ib <- pt[basin_id > 0]
pt_ib[, Basin := sprintf("Basin %d", basin_id)]
basin_order <- sort(unique(pt_ib$Basin))
pt_ib[, Basin := factor(Basin, levels = basin_order)]

fit_ib <- survfit(Surv(time_years, event_death) ~ Basin, data = pt_ib)

n_basins <- length(basin_order)
basin_colors <- brewer.pal(max(3, n_basins), "Set1")[seq_len(n_basins)]
names(basin_colors) <- basin_order

survfit_to_dt <- function(fit) {
  strata_names <- names(fit$strata)
  strata_lens <- fit$strata
  data.table(time = fit$time, surv = fit$surv,
             upper = fit$upper, lower = fit$lower,
             n.risk = fit$n.risk,
             strata = rep(strata_names, strata_lens))
}

sdf <- survfit_to_dt(fit_ib)
sdf[, strata := gsub("Basin=", "", strata)]
time0 <- data.table(time = 0, surv = 1, upper = 1, lower = 1,
                    n.risk = NA_real_, strata = unique(sdf$strata))
sdf <- rbind(time0, sdf)

sd <- survdiff(Surv(time_years, event_death) ~ Basin, data = pt_ib)
pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
pv_label <- if (pv < 0.0001) "Log-rank\np < 0.0001" else sprintf("Log-rank\np = %.4f", pv)

ksz <- 14
g_km <- ggplot(sdf, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12,
              stat = "identity", outline.type = "upper") +
  scale_color_manual(values = basin_colors, name = "") +
  scale_fill_manual(values = basin_colors, name = "") +
  annotate("text", x = 2.5, y = 0.15, label = pv_label,
           size = ksz * 3 / .pt, fontface = "bold", hjust = 0) +
  labs(title = sprintf("%s - Survival by basin", disease_label),
       x = "Years", y = "Survival probability") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 8)) +
  theme_bw(base_size = ksz) +
  theme(plot.title = element_text(size = ksz * 3, face = "bold"),
        axis.text = element_text(size = ksz * 3),
        axis.title = element_text(size = ksz * 3),
        legend.text = element_text(size = ksz * 2.5),
        legend.position = "top",
        legend.key.width = unit(2, "cm"),
        plot.background = element_rect(fill = "white", colour = NA))

risk_times <- seq(0, 8, by = 2)
risk_dt <- data.table()
for (b in basin_order) {
  sub <- pt_ib[Basin == b]
  for (tt in risk_times) {
    risk_dt <- rbind(risk_dt, data.table(
      time = tt, Basin = b, n_risk = sum(sub$time_years >= tt)))
  }
}
risk_dt[, Basin := factor(Basin, levels = rev(basin_order))]
risk_dt[, label := format(n_risk, big.mark = ",")]

g_risk <- ggplot(risk_dt, aes(x = time, y = Basin, label = label, color = Basin)) +
  geom_text(size = ksz * 2.5 / .pt, fontface = "bold") +
  scale_color_manual(values = basin_colors, guide = "none") +
  scale_x_continuous(limits = c(-0.5, 8), breaks = risk_times) +
  labs(x = NULL, y = NULL, title = "Number at risk") +
  theme_bw(base_size = ksz) +
  theme(plot.title = element_text(size = ksz * 2, face = "bold"),
        axis.text.x = element_text(size = ksz * 3),
        axis.text.y = element_text(size = ksz * 2),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 10, 10, 80))

out_km <- file.path(fig_dir, sprintf("Fig%s_b_%s_km.png", fig_num, disease))
png(out_km, width = 24, height = 18, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.75, 0.25), "npc"))))
print(g_km, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g_risk, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
cat(sprintf("  Written: %s\n", out_km))

## ── Panel d: Sankey ──
traj_path <- file.path(base, sprintf("plot/%s/trajectories_dead.tsv", disease))
if (file.exists(traj_path)) {
  traj <- fread(traj_path)
  max_dist <- min(19L, max(traj$periods_before_death))
  hm_data <- traj[periods_before_death <= max_dist]

  hm_data[, group := ifelse(basin_id == 0, "Non-basin", sprintf("Basin %d", basin_id))]

  all_basins_in_data <- sort(unique(hm_data[basin_id > 0, basin_id]))
  all_basin_labels <- sprintf("Basin %d", all_basins_in_data)
  has_nonbasin <- any(hm_data$basin_id == 0)

  sankey_colors <- c(brewer.pal(max(3, length(all_basins_in_data)), "Set1")[seq_along(all_basins_in_data)])
  names(sankey_colors) <- all_basin_labels
  sankey_colors["Non-basin"] <- "gray50"

  legend_breaks <- c(all_basin_labels, "Non-basin")

  all_groups <- sort(unique(hm_data$group))

  sankey_times <- sort(unique(c(seq(max_dist, 0, by = -5), 2, 1, 0)), decreasing = TRUE)

  node_dt <- data.table()
  for (si in seq_along(sankey_times)) {
    t_val <- sankey_times[si]
    at_t <- hm_data[periods_before_death == t_val, .(n = .N), by = group]
    at_t[, total := sum(n)]
    at_t[, frac := n / total]
    for (g in all_groups) {
      if (!g %in% at_t$group)
        at_t <- rbind(at_t, data.table(group = g, n = 0L, total = at_t$total[1], frac = 0))
    }
    setorder(at_t, group)
    y_bottom <- 0
    for (bi in seq_len(nrow(at_t))) {
      node_dt <- rbind(node_dt, data.table(
        x = si, t = t_val, group = at_t$group[bi],
        y_bottom = y_bottom, y_top = y_bottom + at_t$frac[bi],
        n = at_t$n[bi], frac = at_t$frac[bi]))
      y_bottom <- y_bottom + at_t$frac[bi]
    }
  }

  setorder(hm_data, kojin_id, -periods_before_death)
  flow_dt <- data.table()
  for (si in seq_len(length(sankey_times) - 1)) {
    t_from <- sankey_times[si]; t_to <- sankey_times[si + 1]
    pts_from <- hm_data[periods_before_death == t_from, .(kojin_id, from = group)]
    pts_to <- hm_data[periods_before_death == t_to, .(kojin_id, to = group)]
    fl <- merge(pts_from, pts_to, by = "kojin_id")[, .N, by = .(from, to)]
    fl[, total := sum(N)]; fl[, frac := N / total]
    flow_dt <- rbind(flow_dt, data.table(
      x_from = si, x_to = si + 1, t_from = t_from, t_to = t_to,
      from = fl$from, to = fl$to, n = fl$N, frac = fl$frac))
  }

  for (si in seq_len(length(sankey_times) - 1)) {
    fl <- flow_dt[x_from == si]
    if (nrow(fl) == 0) next
    for (b_from in unique(fl$from)) {
      node_from <- node_dt[x == si & group == b_from]
      if (nrow(node_from) == 0 || node_from$frac == 0) next
      sub_fl <- fl[from == b_from]; sub_fl[, flow_frac := n / sum(n)]
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
      sub_fl <- fl[to == b_to]; sub_fl[, flow_frac := n / sum(n)]
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
    x_bot <- rev(x_top); y_bot <- yb1 + (yb2 - yb1) * rev(s)
    data.table(x = c(x_top, x_bot), y = c(y_top, y_bot))
  }

  flow_polys <- list()
  for (i in seq_len(nrow(flow_dt))) {
    r <- flow_dt[i]
    if (is.na(r$y_from_bottom) || r$y_from_top - r$y_from_bottom < 1e-6) next
    poly <- make_band(r$x_from + 0.15, r$x_to - 0.15,
                       r$y_from_bottom, r$y_from_top, r$y_to_bottom, r$y_to_top)
    poly[, `:=`(id = i, from = r$from, to = r$to)]
    flow_polys[[i]] <- poly
  }
  flow_poly_dt <- rbindlist(flow_polys)
  flow_poly_dt[, color := ifelse(from == to, from, "transition")]

  time_labels <- ifelse(sankey_times == 0, "death", as.character(sankey_times))

  dummy_rect <- data.table(xmin = -99, xmax = -98, ymin = -99, ymax = -98, group = "Non-basin")

  g_s <- ggplot() +
    geom_rect(data = dummy_rect,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group)) +
    geom_polygon(data = flow_poly_dt,
                 aes(x = x, y = y, group = id, fill = color), alpha = 0.4) +
    geom_rect(data = node_dt[frac > 0],
              aes(xmin = x - 0.12, xmax = x + 0.12,
                  ymin = y_bottom, ymax = y_top, fill = group),
              color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c(sankey_colors, "transition" = "gray70"),
                      name = NULL, breaks = legend_breaks, drop = FALSE) +
    coord_cartesian(xlim = c(0.5, length(sankey_times) + 0.5), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq_along(sankey_times), labels = time_labels) +
    scale_y_continuous(labels = NULL) +
    labs(title = sprintf("%s: basin transitions before death (n=%s)",
                          disease_label, format(uniqueN(hm_data$kojin_id), big.mark = ",")),
         x = "Periods before death", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 33),
          axis.title.x = element_text(size = 33))

  out_s <- file.path(fig_dir, sprintf("Fig%s_d_%s_sankey.png", fig_num, disease))
  ggsave(out_s, g_s, width = 14, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  Written: %s\n", out_s))
}

## ── Panel c: copy existing logitdiff barcharts ──
for (bid in all_basins_in_data) {
  src <- file.path(base, sprintf("plot/%s/tagcloud_scored/logitdiff/barchart_logitdiff_basin%d.png", disease, bid))
  if (file.exists(src)) {
    dst <- file.path(fig_dir, sprintf("barchart_logitdiff_basin%d.png", bid))
    file.copy(src, dst, overwrite = TRUE)
    cat(sprintf("  Copied: %s\n", dst))
  }
}

cat(sprintf("\n=== Fig%s (%s) complete ===\n", fig_num, disease_label))
