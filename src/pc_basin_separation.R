## Systematic PC-basin separation analysis across all disease families
##
## For each disease x PC: compute basin separation score (eta-squared from ANOVA)
## Identify the most basin-discriminating PC per disease, extract top loadings
##
## Usage:
##   Rscript src/pc_basin_separation.R <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

all_scores <- list()
all_loadings <- list()

for (fid in fids) {
  scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  eigvec_path <- sprintf("output/%s/all_comorbid/pca/7/Eigen_vectors.csv", fid)
  col_ids_path <- sprintf("output/%s/all_comorbid/cohort/col_ids.txt", fid)

  if (!file.exists(scores_path) || !file.exists(eigvec_path)) next

  pt <- fread(scores_path)
  if (is.null(pt$basin_id) || length(unique(pt$basin_id)) < 2) next

  pt <- pt[!is.na(basin_id)]
  pt[, basin_f := factor(basin_id)]
  n_basins <- length(unique(pt$basin_id))

  # --- Eta-squared for each PC ---
  pc_cols <- paste0("PC", 1:7)
  pc_cols <- intersect(pc_cols, names(pt))

  for (pc in pc_cols) {
    vals <- pt[[pc]]
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) next
    ss_total <- sum((vals - mean(vals, na.rm = TRUE))^2, na.rm = TRUE)
    group_means <- pt[, .(m = mean(get(pc), na.rm = TRUE), n = .N), by = basin_f]
    ss_between <- sum(group_means$n * (group_means$m - mean(vals, na.rm = TRUE))^2)
    eta2 <- ss_between / ss_total

    # Correlation with energy
    cor_energy <- cor(vals, pt$energy, use = "pairwise.complete.obs")

    # Correlation with death
    pt[, death_bin := as.integer(death == "dead")]
    cor_death <- cor(vals, pt$death_bin, use = "pairwise.complete.obs")

    all_scores[[length(all_scores) + 1]] <- data.table(
      disease = fid, pc = pc, n_basins = n_basins,
      eta2_basin = eta2,
      cor_energy = cor_energy,
      cor_death = cor_death,
      n_patients = nrow(pt)
    )
  }

  # --- Top loadings for each PC ---
  if (!file.exists(col_ids_path)) next
  eigvec <- tryCatch(as.matrix(read.csv(eigvec_path, header = FALSE)),
                     error = function(e) NULL)
  if (is.null(eigvec)) next

  col_lines <- readLines(col_ids_path, encoding = "UTF-8")
  code_names <- character(nrow(eigvec))
  for (ln in col_lines) {
    parts <- strsplit(ln, ",")[[1]]
    if (length(parts) >= 2) {
      code <- trimws(parts[1])
      id <- as.integer(trimws(parts[2]))
      if (!is.na(id) && id >= 1 && id <= length(code_names)) {
        code_names[id] <- code
      }
    }
  }

  n_top <- 10
  for (pc_idx in seq_along(pc_cols)) {
    if (pc_idx > ncol(eigvec)) break
    v <- eigvec[, pc_idx]
    top_idx <- order(-abs(v))[1:min(n_top, length(v))]
    for (ti in top_idx) {
      all_loadings[[length(all_loadings) + 1]] <- data.table(
        disease = fid, pc = pc_cols[pc_idx],
        rank = which(top_idx == ti),
        icd10 = code_names[ti],
        loading = v[ti]
      )
    }
  }

  cat(sprintf("  %s: %d basins, done\n", fid, n_basins))
}

scores_dt <- rbindlist(all_scores)
loadings_dt <- rbindlist(all_loadings)

# --- Best PC per disease (highest eta2) ---
best_pc <- scores_dt[, .SD[which.max(eta2_basin)], by = disease]
best_pc <- best_pc[order(-eta2_basin)]

fwrite(scores_dt, file.path(out_dir, "pc_basin_scores_all.tsv"), sep = "\t")
fwrite(best_pc, file.path(out_dir, "pc_basin_best.tsv"), sep = "\t")
fwrite(loadings_dt, file.path(out_dir, "pc_top_loadings_all.tsv"), sep = "\t")

cat("\n=== Disease ranking by best PC basin separation (eta-squared) ===\n")
print(best_pc[, .(disease, pc, eta2_basin, cor_energy, cor_death, n_basins)], nrows = 60)

# --- Plot: eta-squared ranking ---
best_pc[, disease := factor(disease, levels = rev(disease))]
g1 <- ggplot(best_pc, aes(x = eta2_basin, y = disease, fill = pc)) +
  geom_col() +
  labs(title = "Basin separation by best PC (eta-squared)",
       subtitle = "Higher = PC more strongly distinguishes basins",
       x = expression(eta^2), y = NULL, fill = "Best PC") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "pc_basin_separation_ranking.png"),
       g1, width = 10, height = 14, dpi = 150)

# --- Plot: all PCs heatmap per disease ---
scores_dt[, disease := factor(disease, levels = rev(best_pc$disease))]
g2 <- ggplot(scores_dt, aes(x = pc, y = disease, fill = eta2_basin)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#1565C0") +
  labs(title = "PC x Disease: Basin separation (eta-squared)",
       x = NULL, y = NULL, fill = expression(eta^2)) +
  theme_minimal(base_size = 9)
ggsave(file.path(out_dir, "pc_basin_heatmap.png"),
       g2, width = 8, height = 14, dpi = 150)

# --- Extract loadings for best PC per disease ---
best_loadings <- merge(best_pc[, .(disease, pc)], loadings_dt,
                       by = c("disease", "pc"))
fwrite(best_loadings, file.path(out_dir, "best_pc_loadings.tsv"), sep = "\t")

cat(sprintf("\nOutput: %s\n", out_dir))
cat("  pc_basin_scores_all.tsv\n")
cat("  pc_basin_best.tsv\n")
cat("  pc_top_loadings_all.tsv\n")
cat("  best_pc_loadings.tsv\n")
cat("  pc_basin_separation_ranking.png\n")
cat("  pc_basin_heatmap.png\n")
