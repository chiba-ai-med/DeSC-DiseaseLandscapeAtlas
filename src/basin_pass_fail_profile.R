## Compare PASS vs FAIL basins: death rate, age, energy, sex
## Show that FAIL basins are "milder" (lower death rate, younger, etc.)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "summary"

# Collect basin-stratified Cox results
cox_files <- Sys.glob("plot/*/cox_basin_stratified.tsv")
cox_all <- rbindlist(lapply(cox_files, fread))

# Identify diseases with both PASS and FAIL basins
mixed <- cox_all[, .(has_pass = any(direction_ok), has_fail = any(!direction_ok)),
                 by = disease_family]
mixed <- mixed[has_pass == TRUE & has_fail == TRUE, disease_family]
cat(sprintf("Diseases with mixed PASS/FAIL basins: %d\n", length(mixed)))
cat(paste(mixed, collapse = ", "), "\n\n")

# For each mixed disease, compute basin profiles
rows <- list()

for (fid in sort(unique(cox_all$disease_family))) {
  scores_path <- sprintf("output/%s/all_comorbid/patient/patient_scores.tsv", fid)
  if (!file.exists(scores_path)) next

  pt <- fread(scores_path)
  pt <- pt[!is.na(basin_id) & is.finite(energy)]
  pt[, death_bin := as.integer(death == "dead")]
  pt[, sex_f := as.integer(sex == "female")]

  cox_sub <- cox_all[disease_family == fid]

  for (bid in unique(pt$basin_id)) {
    bp <- pt[basin_id == bid]
    cox_row <- cox_sub[basin_id == bid]
    if (nrow(cox_row) == 0) next

    rows[[length(rows) + 1]] <- data.table(
      disease = fid,
      basin_id = bid,
      n_patients = nrow(bp),
      death_rate = mean(bp$death_bin, na.rm = TRUE),
      mean_age = mean(bp$age, na.rm = TRUE),
      mean_energy = mean(bp$energy, na.rm = TRUE),
      sd_energy = sd(bp$energy, na.rm = TRUE),
      pct_female = mean(bp$sex_f, na.rm = TRUE),
      HR = cox_row$HR,
      pass = cox_row$direction_ok
    )
  }
}

dt <- rbindlist(rows)
dt[, status := ifelse(pass, "PASS", "FAIL")]

fwrite(dt, file.path(out_dir, "basin_profiles.tsv"), sep = "\t")

# --- Mixed diseases only ---
dt_mixed <- dt[disease %in% mixed]
cat("=== Mixed PASS/FAIL basin profiles ===\n")
print(dt_mixed[order(disease, basin_id),
               .(disease, basin_id, status, n_patients, death_rate,
                 mean_age, mean_energy, pct_female)], nrows = 100)

# --- Aggregate: PASS vs FAIL basins ---
agg <- dt[, .(
  n_basins = .N,
  median_death_rate = median(death_rate),
  mean_death_rate = mean(death_rate),
  median_age = median(mean_age),
  median_energy = median(mean_energy),
  median_pct_female = median(pct_female),
  n_diseases = length(unique(disease))
), by = status]
cat("\n=== Aggregate: PASS vs FAIL basins ===\n")
print(agg)

# --- Test: FAIL basins have lower death rate? ---
pass_dr <- dt[pass == TRUE, death_rate]
fail_dr <- dt[pass == FALSE, death_rate]
wt <- wilcox.test(pass_dr, fail_dr)
cat(sprintf("\nDeath rate: PASS median=%.4f, FAIL median=%.4f, Wilcoxon p=%.2e\n",
            median(pass_dr), median(fail_dr), wt$p.value))

pass_age <- dt[pass == TRUE, mean_age]
fail_age <- dt[pass == FALSE, mean_age]
wt2 <- wilcox.test(pass_age, fail_age)
cat(sprintf("Age: PASS median=%.1f, FAIL median=%.1f, Wilcoxon p=%.2e\n",
            median(pass_age), median(fail_age), wt2$p.value))

pass_e <- dt[pass == TRUE, mean_energy]
fail_e <- dt[pass == FALSE, mean_energy]
wt3 <- wilcox.test(pass_e, fail_e)
cat(sprintf("Energy: PASS median=%.3f, FAIL median=%.3f, Wilcoxon p=%.2e\n",
            median(pass_e), median(fail_e), wt3$p.value))

# --- Plot: death rate PASS vs FAIL ---
dt_mixed[, label := paste0(disease, " B", basin_id)]
dt_mixed[, label := factor(label, levels = rev(label[order(death_rate)]))]

g1 <- ggplot(dt_mixed, aes(x = death_rate, y = label, fill = status)) +
  geom_col() +
  scale_fill_manual(values = c(PASS = "#1565C0", FAIL = "#C62828")) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Death rate by basin (mixed PASS/FAIL diseases)",
       subtitle = "FAIL basins tend to have lower death rates",
       x = "Death rate", y = NULL, fill = NULL) +
  theme_minimal(base_size = 9)
ggsave(file.path(out_dir, "basin_pass_fail_death_rate.png"),
       g1, width = 10, height = 12, dpi = 150)

# --- Plot: paired comparison within each disease ---
dt_mixed_wide <- dcast(dt_mixed, disease ~ status,
                       value.var = "death_rate", fun.aggregate = mean)
dt_mixed_wide <- dt_mixed_wide[!is.na(PASS) & !is.na(FAIL)]

g2 <- ggplot(dt_mixed_wide, aes(x = PASS, y = FAIL)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = disease), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray40") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Death rate: PASS basin vs FAIL basin (within same disease)",
       subtitle = "Points below diagonal = FAIL basin has lower death rate",
       x = "PASS basin death rate", y = "FAIL basin death rate") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "basin_pass_vs_fail_scatter.png"),
       g2, width = 9, height = 8, dpi = 150)

cat("\nDone.\n")
