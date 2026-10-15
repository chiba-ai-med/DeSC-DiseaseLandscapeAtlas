## PC axis correlations with age, sex, energy for all disease families
##
## For each disease x PC: cor(PC, age), cor(PC, sex_f), cor(PC, energy)
## Output: correlation table + heatmap

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "summary"
fids <- list.dirs("output", recursive = FALSE, full.names = FALSE)

rows <- list()

for (fid in fids) {
  scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  if (!file.exists(scores_path)) next

  pt <- fread(scores_path)
  pt <- pt[is.finite(energy)]
  if (nrow(pt) < 100) next

  pt[, sex_num := as.integer(sex == "female")]

  pc_cols <- paste0("PC", 1:7)
  pc_cols <- intersect(pc_cols, names(pt))

  for (pc in pc_cols) {
    vals <- pt[[pc]]
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) next

    rows[[length(rows) + 1]] <- data.table(
      disease = fid, pc = pc,
      cor_age = cor(vals, pt$age, use = "pairwise.complete.obs"),
      cor_sex = cor(vals, pt$sex_num, use = "pairwise.complete.obs"),
      cor_energy = cor(vals, pt$energy, use = "pairwise.complete.obs")
    )
  }

  # energy vs age, energy vs sex
  rows[[length(rows) + 1]] <- data.table(
    disease = fid, pc = "Energy",
    cor_age = cor(pt$energy, pt$age, use = "pairwise.complete.obs"),
    cor_sex = cor(pt$energy, pt$sex_num, use = "pairwise.complete.obs"),
    cor_energy = 1.0
  )

  cat(sprintf("  %s\n", fid))
}

dt <- rbindlist(rows)
fwrite(dt, file.path(out_dir, "pc_age_sex_energy_correlations.tsv"), sep = "\t")

# --- Energy vs age summary ---
energy_dt <- dt[pc == "Energy"]
energy_dt <- energy_dt[order(cor_age)]
cat("\n=== Energy-Age correlation by disease ===\n")
print(energy_dt[, .(disease, cor_age, cor_sex)], nrows = 60)

cat(sprintf("\nEnergy-Age: median=%.3f, range=[%.3f, %.3f]\n",
            median(energy_dt$cor_age), min(energy_dt$cor_age), max(energy_dt$cor_age)))
cat(sprintf("Energy-Sex: median=%.3f, range=[%.3f, %.3f]\n",
            median(energy_dt$cor_sex), min(energy_dt$cor_sex), max(energy_dt$cor_sex)))

# --- Heatmap: disease x PC, colored by cor_age ---
pc_dt <- dt[pc != "Energy"]
disease_order <- energy_dt[order(cor_age), disease]
pc_dt[, disease_f := factor(disease, levels = rev(disease_order))]

g1 <- ggplot(pc_dt, aes(x = pc, y = disease_f, fill = cor_age)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#C62828", mid = "white", high = "#1565C0",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(title = "cor(PC, Age) by disease family",
       x = NULL, y = NULL, fill = "cor") +
  theme_minimal(base_size = 9)

g2 <- ggplot(pc_dt, aes(x = pc, y = disease_f, fill = cor_sex)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#C62828", mid = "white", high = "#1565C0",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(title = "cor(PC, Sex[female]) by disease family",
       x = NULL, y = NULL, fill = "cor") +
  theme_minimal(base_size = 9)

g3 <- ggplot(pc_dt, aes(x = pc, y = disease_f, fill = cor_energy)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#C62828", mid = "white", high = "#1565C0",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(title = "cor(PC, Energy) by disease family",
       x = NULL, y = NULL, fill = "cor") +
  theme_minimal(base_size = 9)

ggsave(file.path(out_dir, "pc_cor_age_heatmap.png"), g1, width = 8, height = 14, dpi = 150)
ggsave(file.path(out_dir, "pc_cor_sex_heatmap.png"), g2, width = 8, height = 14, dpi = 150)
ggsave(file.path(out_dir, "pc_cor_energy_heatmap.png"), g3, width = 8, height = 14, dpi = 150)

# --- Energy vs age/sex scatter across diseases ---
g4 <- ggplot(energy_dt, aes(x = cor_age, y = cor_sex)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = disease), size = 2.5, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Energy correlation: Age vs Sex (by disease family)",
       x = "cor(Energy, Age)", y = "cor(Energy, Sex[female])") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "energy_age_sex_scatter.png"), g4, width = 10, height = 8, dpi = 150)

cat("\nDone.\n")
