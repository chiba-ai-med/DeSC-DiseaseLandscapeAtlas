## PCA scatter plots (PC1-7, all diseases combined) colored by attributes
##
## Usage:
##   Rscript src/plot_pca_pairs_all.R <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load all diseases ---
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
cat(sprintf("Total: %d patients from %d diseases\n", nrow(dt), length(all_pts)))

dt[, death_bin := as.integer(death == "dead")]
dt[, sex_f := as.integer(sex == "female")]
dt[, kenshin := as.integer(kenshin_data_ari)]
dt[, hosp_bin := as.integer(hospitalized == "yes")]
dt[is.na(hosp_bin), hosp_bin := as.integer(hospitalized)]

set.seed(42)
n_sample <- min(200000L, nrow(dt))
pt_sub <- dt[sample.int(nrow(dt), n_sample)]
rm(dt)
cat(sprintf("Sampled: %d patients\n", n_sample))

# --- PC pairs (PC1-7: 21 pairs) ---
pc_cols <- paste0("PC", 1:7)
pc_pairs <- list()
for (i in 1:6) {
  for (j in (i + 1):7) {
    pc_pairs[[length(pc_pairs) + 1]] <- c(pc_cols[i], pc_cols[j])
  }
}
cat(sprintf("PC pairs: %d\n", length(pc_pairs)))

# --- Color variables ---
labels <- list(
  age          = list(col = "age",          type = "cont", title = "Age"),
  sex          = list(col = "sex_f",        type = "cont", title = "Sex (Female=1)"),
  death        = list(col = "death_bin",    type = "cont", title = "Death"),
  energy       = list(col = "energy",       type = "cont", title = "Energy"),
  hospitalized = list(col = "hosp_bin",     type = "cont", title = "Hospitalized"),
  n_admissions = list(col = "n_admissions", type = "cont", title = "N admissions"),
  basin        = list(col = "basin_id",     type = "cat",  title = "Basin"),
  kenshin      = list(col = "kenshin",      type = "cont", title = "Kenshin"),
  obs_years    = list(col = "obs_years",    type = "cont", title = "Obs years"),
  n_diagnoses  = list(col = "n_diagnoses",  type = "cont", title = "N diagnoses")
)

for (lname in names(labels)) {
  ldef <- labels[[lname]]
  col <- ldef$col
  if (!(col %in% names(pt_sub))) next

  sub <- pt_sub[is.finite(get(col))]
  if (nrow(sub) < 1000) {
    cat(sprintf("  Skipping %s (too few valid rows: %d)\n", lname, nrow(sub)))
    next
  }
  cat(sprintf("  Plotting %s (n=%d)...\n", lname, nrow(sub)))

  plots <- list()
  for (pi in seq_along(pc_pairs)) {
    pcx <- pc_pairs[[pi]][1]
    pcy <- pc_pairs[[pi]][2]
    d <- sub[is.finite(get(pcx)) & is.finite(get(pcy))]

    if (ldef$type == "cat") {
      d[[col]] <- factor(d[[col]])
      p <- ggplot(d, aes(x = .data[[pcx]], y = .data[[pcy]],
                          color = .data[[col]])) +
        geom_point(size = 0.08, alpha = 0.15) +
        labs(x = pcx, y = pcy, color = ldef$title) +
        theme_minimal(base_size = 7) +
        theme(legend.position = if (pi == 1) "right" else "none",
              plot.margin = margin(1, 1, 1, 1))
    } else {
      p <- ggplot(d, aes(x = .data[[pcx]], y = .data[[pcy]],
                          color = .data[[col]])) +
        geom_point(size = 0.08, alpha = 0.15) +
        scale_color_viridis_c(option = "C") +
        labs(x = pcx, y = pcy, color = ldef$title) +
        theme_minimal(base_size = 7) +
        theme(legend.position = if (pi == 1) "right" else "none",
              plot.margin = margin(1, 1, 1, 1))
    }
    plots[[pi]] <- p
  }

  outfile <- file.path(out_dir, sprintf("pairs_all_%s.png", lname))
  n_cols <- 7
  n_rows <- ceiling(length(plots) / n_cols)

  png(outfile, width = 300 * n_cols, height = 280 * n_rows + 40, res = 150)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows, n_cols)))
  for (pi in seq_along(plots)) {
    row <- ((pi - 1) %/% n_cols) + 1
    col_pos <- ((pi - 1) %% n_cols) + 1
    print(plots[[pi]],
          vp = grid::viewport(layout.pos.row = row, layout.pos.col = col_pos))
  }
  grid::grid.text(
    sprintf("All diseases: PC1-7 pairs by %s (n=%d)", ldef$title, nrow(sub)),
    y = 0.99, gp = grid::gpar(fontsize = 10))
  dev.off()
  cat(sprintf("    -> %s\n", outfile))
}

cat("Done\n")
