## PC pair plots colored by attributes (all diseases combined)
## Hex bins with mean attribute value as color

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "summary/pc_pairplots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

all_pts <- list()
for (fid in fids) {
  path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  if (!file.exists(path)) next
  pt <- fread(path)
  pt[, disease := fid]
  ndiag_path <- sprintf("output/%s/all_comorbid/cohort/n_diagnoses.tsv", fid)
  if (file.exists(ndiag_path)) {
    nd <- fread(ndiag_path)
    pt <- merge(pt, nd, by = "kojin_id", all.x = TRUE)
  } else {
    pt[, n_diagnoses := NA_integer_]
  }
  all_pts[[length(all_pts) + 1]] <- pt
}
dt <- rbindlist(all_pts, fill = TRUE)
dt[, death_bin := as.integer(death == "dead")]
dt[, sex_f := as.integer(sex == "female")]
dt[, kenshin := as.integer(kenshin_data_ari)]
dt[, hosp_bin := as.integer(hospitalized == "yes")]
dt[is.na(hosp_bin), hosp_bin := as.integer(hospitalized)]

cat(sprintf("Total: %d patients\n", nrow(dt)))

set.seed(42)
samp <- dt[sample(.N, min(200000, .N))]

pc_cols <- paste0("PC", 1:7)
pc_pairs <- list()
for (i in 1:6) {
  for (j in (i + 1):7) {
    pc_pairs[[length(pc_pairs) + 1]] <- c(pc_cols[i], pc_cols[j])
  }
}

color_vars <- list(
  list(name = "energy", col = "energy", label = "Energy", palette = "RdBu"),
  list(name = "age", col = "age", label = "Age", palette = "viridis"),
  list(name = "sex", col = "sex_f", label = "Sex (Female)", palette = "RdBu"),
  list(name = "death", col = "death_bin", label = "Death", palette = "Reds"),
  list(name = "kenshin", col = "kenshin", label = "Kenshin", palette = "Blues"),
  list(name = "hospitalized", col = "hosp_bin", label = "Hospitalized", palette = "Reds"),
  list(name = "n_admissions", col = "n_admissions", label = "N admissions", palette = "YlOrRd"),
  list(name = "obs_years", col = "obs_years", label = "Obs years", palette = "viridis")
)

for (cv in color_vars) {
  cat(sprintf("  Plotting %s...\n", cv$name))
  plots <- list()
  for (pair in pc_pairs) {
    pcx <- pair[1]
    pcy <- pair[2]
    sub <- samp[is.finite(get(cv$col))]

    g <- ggplot(sub, aes(x = .data[[pcx]], y = .data[[pcy]],
                          z = .data[[cv$col]])) +
      stat_summary_hex(fun = mean, bins = 60) +
      labs(x = pcx, y = pcy, fill = cv$label) +
      theme_minimal(base_size = 7) +
      theme(legend.position = "none",
            plot.margin = margin(2, 2, 2, 2))

    if (cv$palette == "RdBu") {
      g <- g + scale_fill_gradient2(low = "#C62828", mid = "white",
                                     high = "#1565C0", midpoint = mean(sub[[cv$col]], na.rm = TRUE))
    } else if (cv$palette == "viridis") {
      g <- g + scale_fill_viridis_c()
    } else if (cv$palette == "Reds") {
      g <- g + scale_fill_gradient(low = "white", high = "#C62828")
    } else if (cv$palette == "Blues") {
      g <- g + scale_fill_gradient(low = "white", high = "#1565C0")
    } else if (cv$palette == "YlOrRd") {
      g <- g + scale_fill_gradient(low = "#FFFFCC", high = "#C62828")
    }
    plots[[length(plots) + 1]] <- g
  }

  # Arrange in a triangular grid: 6 rows, up to 6 cols
  # pair (i,j): row = i, col = j-i
  n_pc <- 7
  blank <- ggplot() + theme_void()

  grid <- list()
  idx <- 1
  for (i in 1:(n_pc - 1)) {
    for (j in 1:(n_pc - 1)) {
      if (j > i) {
        grid[[length(grid) + 1]] <- blank
      } else if (j == i) {
        # Diagonal: just label
        grid[[length(grid) + 1]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = pc_cols[i + 1],
                   size = 4, fontface = "bold") +
          theme_void()
      } else {
        grid[[length(grid) + 1]] <- plots[[idx]]
        idx <- idx + 1
      }
    }
  }

  # Actually simpler: just arrange the 21 plots in a grid
  ncols <- 7
  nrows <- 3
  out_file <- file.path(out_dir, sprintf("pc_pairs_%s.png", cv$name))

  # Use a simple cowplot-free approach: save each plot and note it
  # Or: just arrange 21 plots in 3x7
  pdf(NULL)  # prevent empty Rplots.pdf
  png(out_file, width = 2100, height = 900, res = 150)
  par(mfrow = c(1, 1))
  dev.off()

  # Better: use ggplot faceting with a long-format approach
  pair_labels <- sapply(pc_pairs, function(p) paste(p[1], "vs", p[2]))
  long_list <- list()
  for (k in seq_along(pc_pairs)) {
    pair <- pc_pairs[[k]]
    tmp <- samp[is.finite(get(cv$col)), .(x = get(pair[1]), y = get(pair[2]),
                                            z = get(cv$col))]
    tmp[, panel := pair_labels[k]]
    long_list[[k]] <- tmp
  }
  long_dt <- rbindlist(long_list)
  long_dt[, panel := factor(panel, levels = pair_labels)]

  g_all <- ggplot(long_dt, aes(x = x, y = y, z = z)) +
    stat_summary_hex(fun = mean, bins = 40) +
    facet_wrap(~ panel, scales = "free", ncol = 7) +
    labs(title = sprintf("PC pair plots colored by %s (all diseases, n=200k sample)",
                         cv$label),
         x = NULL, y = NULL, fill = cv$label) +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(size = 6))

  if (cv$palette == "RdBu") {
    g_all <- g_all + scale_fill_gradient2(low = "#C62828", mid = "white",
                                           high = "#1565C0",
                                           midpoint = mean(long_dt$z, na.rm = TRUE))
  } else if (cv$palette == "viridis") {
    g_all <- g_all + scale_fill_viridis_c()
  } else if (cv$palette == "Reds") {
    g_all <- g_all + scale_fill_gradient(low = "white", high = "#C62828")
  } else if (cv$palette == "Blues") {
    g_all <- g_all + scale_fill_gradient(low = "white", high = "#1565C0")
  } else if (cv$palette == "YlOrRd") {
    g_all <- g_all + scale_fill_gradient(low = "#FFFFCC", high = "#C62828")
  }

  ggsave(out_file, g_all, width = 21, height = 9, dpi = 150)
}

cat("\nDone.\n")
