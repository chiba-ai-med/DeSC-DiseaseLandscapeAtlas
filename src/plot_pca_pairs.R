## PCA scatter plots with label overlay (ggplot2 only)
##
## Usage:
##   Rscript src/plot_pca_pairs.R <patient_scores_tsv> <disease> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
scores_path <- args[1]
disease     <- args[2]
out_dir     <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("=== plot_pca_pairs: %s ===\n", disease))

pt <- fread(scores_path)
cat(sprintf("  %d patients\n", nrow(pt)))

set.seed(42)
n_sample <- min(30000L, nrow(pt))
pt_sub <- pt[sample.int(nrow(pt), n_sample)]
rm(pt)

labels <- list(
  age          = list(col = "age",          type = "cont", title = "Age"),
  sex          = list(col = "sex",          type = "cat",  title = "Sex"),
  death        = list(col = "death",        type = "cat",  title = "Death"),
  hospitalized = list(col = "hospitalized", type = "cat",  title = "Hospitalization"),
  n_admissions = list(col = "n_admissions", type = "cont", title = "N admissions"),
  basin        = list(col = "basin",        type = "cat",  title = "Basin"),
  energy       = list(col = "energy",       type = "cont", title = "Energy"),
  insurer      = list(col = "insurer",      type = "cat",  title = "Insurer")
)
labels <- labels[vapply(labels, function(l) l$col %in% names(pt_sub), logical(1))]
cat(sprintf("  Labels: %s\n", paste(names(labels), collapse = ", ")))

pc_pairs <- list(
  c("PC1", "PC2"), c("PC1", "PC3"), c("PC1", "PC4"),
  c("PC2", "PC3"), c("PC2", "PC4"), c("PC3", "PC4")
)

for (lname in names(labels)) {
  ldef <- labels[[lname]]
  col <- ldef$col
  cat(sprintf("  Plotting %s ...\n", lname))

  plots <- list()
  for (pi in seq_along(pc_pairs)) {
    pcx <- pc_pairs[[pi]][1]
    pcy <- pc_pairs[[pi]][2]
    d <- pt_sub[is.finite(get(pcx)) & is.finite(get(pcy))]

    if (ldef$type == "cat") {
      d[[col]] <- factor(d[[col]])
      p <- ggplot(d, aes(x = .data[[pcx]], y = .data[[pcy]],
                          color = .data[[col]])) +
        geom_point(size = 0.15, alpha = 0.25) +
        labs(x = pcx, y = pcy, color = ldef$title) +
        theme_minimal(base_size = 9) +
        theme(legend.position = if (pi == 1) "right" else "none")
    } else {
      d <- d[is.finite(get(col))]
      p <- ggplot(d, aes(x = .data[[pcx]], y = .data[[pcy]],
                          color = .data[[col]])) +
        geom_point(size = 0.15, alpha = 0.25) +
        scale_color_viridis_c(option = "C") +
        labs(x = pcx, y = pcy, color = ldef$title) +
        theme_minimal(base_size = 9) +
        theme(legend.position = if (pi == 1) "right" else "none")
    }
    plots[[pi]] <- p
  }

  outfile <- file.path(out_dir, sprintf("pairs_%s_%s.png", disease, lname))
  png(outfile, width = 1800, height = 1000, res = 150)
  par(mfrow = c(1, 1))
  # Use viewport layout from grid
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
  for (pi in seq_along(plots)) {
    row <- ((pi - 1) %/% 3) + 1
    col_pos <- ((pi - 1) %% 3) + 1
    print(plots[[pi]],
          vp = grid::viewport(layout.pos.row = row, layout.pos.col = col_pos))
  }
  grid::grid.text(
    sprintf("%s: PC1-4 pairs by %s (n=%d)", disease, ldef$title, n_sample),
    y = 0.99, gp = grid::gpar(fontsize = 13))
  dev.off()
  cat(sprintf("    -> %s\n", outfile))
}

cat("Done\n")
