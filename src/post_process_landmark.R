## Post-process landmark Cox results:
## 1. Compute per-basin majority criterion
## 2. Update atlas_summary.tsv with per-basin majority direction_ok
## 3. Regenerate atlas heatmap

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)
SCRATCH <- "/tmp/claude-1002/-home-koki-dev-DeSC-DiseaseLandscapeAtlas-paper-claude/f5e2f2fe-00ad-43ab-ac3c-a2358eb43a33/scratchpad"

# ============================================================
# 1. Load results
# ============================================================
overall <- fread(file.path(SCRATCH, "landmark_all_diseases_overall.tsv"))
per_basin <- fread(file.path(SCRATCH, "landmark_all_diseases_per_basin.tsv"))

cat(sprintf("Overall: %d diseases\n", nrow(overall)))
cat(sprintf("Per-basin: %d rows across %d diseases\n",
            nrow(per_basin), length(unique(per_basin$disease))))

# ============================================================
# 2. Per-basin majority criterion
# ============================================================
pb_summary <- per_basin[basin_id > 0 & direction %in% c("PASS", "FAIL"),
                        .(n_basins_valid = .N,
                          n_pass = sum(direction == "PASS")),
                        by = disease]
pb_summary[, majority_pass := n_pass >= n_basins_valid / 2]

# Single-basin diseases: use overall direction
single_basin <- overall[!disease %in% pb_summary$disease & direction %in% c("PASS", "FAIL")]
if (nrow(single_basin) > 0) {
  single_basin[, `:=`(n_basins_valid = 1L,
                       n_pass = as.integer(direction == "PASS"),
                       majority_pass = direction == "PASS")]
  pb_summary <- rbind(pb_summary, single_basin[, .(disease, n_basins_valid, n_pass, majority_pass)])
}

cat("\n=== Per-basin majority criterion ===\n")
cat(sprintf("  Diseases with majority PASS: %d\n", sum(pb_summary$majority_pass)))
cat(sprintf("  Diseases with majority FAIL: %d\n", sum(!pb_summary$majority_pass)))

# Detailed per-basin results
cat("\n--- Per-basin detail ---\n")
for (d in sort(pb_summary$disease)) {
  pb_d <- per_basin[disease == d & basin_id > 0 & direction %in% c("PASS", "FAIL")]
  ov_d <- overall[disease == d]
  maj <- pb_summary[disease == d, majority_pass]
  if (nrow(pb_d) > 0) {
    detail <- paste(sprintf("B%s:HR=%.3f(%s)", pb_d$basin_id, pb_d$HR, pb_d$direction), collapse=", ")
  } else {
    detail <- sprintf("overall:HR=%.3f(%s)", ov_d$HR, ov_d$direction)
  }
  cat(sprintf("  %-35s %d/%d → %s  [%s]\n",
              d,
              pb_summary[disease == d, n_pass],
              pb_summary[disease == d, n_basins_valid],
              ifelse(maj, "PASS", "FAIL"),
              detail))
}

# ============================================================
# 3. Update atlas_summary.tsv with per-basin majority criterion
# ============================================================
cat("\n=== Updating atlas_summary.tsv ===\n")

atlas <- fread("summary/atlas_summary.tsv")
cat(sprintf("  Existing atlas rows: %d\n", nrow(atlas)))

# Merge per-basin majority
# Compute per-basin majority for pipeline-compatible output
# The per_basin results include basin_id=0 (non-basin); exclude it
per_basin_valid <- per_basin[basin_id > 0 & direction %in% c("PASS", "FAIL")]
pb_maj <- per_basin_valid[, .(n_pass = sum(direction == "PASS"),
                               n_total = .N),
                           by = disease]
pb_maj[, majority_pass := n_pass >= n_total / 2]

atlas[, direction_ok_overall := direction_ok]
atlas[, direction_ok := NULL]

maj_dt <- pb_maj[, .(disease_family = disease, majority_pass)]
atlas <- merge(atlas, maj_dt, by = "disease_family", all.x = TRUE)
atlas[, direction_ok := majority_pass]
atlas[, majority_pass := NULL]

atlas[is.na(direction_ok), direction_ok := direction_ok_overall]
atlas[, direction_ok_overall := NULL]

# Also update per-disease pipeline output files to exclude basin_id=0
for (d in unique(per_basin_valid$disease)) {
  out_dir <- file.path("output", d, "all_comorbid", "cox", "7")
  if (!dir.exists(out_dir)) next
  pb_d <- per_basin_valid[disease == d]
  if (nrow(pb_d) == 0) next
  pb_out <- data.table(
    disease_family = d,
    basin_id = pb_d$basin_id,
    outcome = "mortality",
    HR = pb_d$HR,
    HR_lower = pb_d$HR_lo,
    HR_upper = pb_d$HR_hi,
    pvalue = pb_d$pval,
    direction_ok = pb_d$direction == "PASS",
    n_patients = pb_d$N,
    n_deaths = pb_d$deaths,
    death_rate = pb_d$deaths / pb_d$N,
    note = ""
  )
  fwrite(pb_out, file.path(out_dir, "energy_sanity_check_per_basin.tsv"), sep = "\t")
}

fwrite(atlas, "summary/atlas_summary.tsv", sep = "\t")
cat(sprintf("  Updated atlas_summary.tsv: %d rows\n", nrow(atlas)))

# Summary counts
n_pass <- sum(atlas$direction_ok == TRUE, na.rm = TRUE)
n_fail <- sum(atlas$direction_ok == FALSE, na.rm = TRUE)
n_na <- sum(is.na(atlas$direction_ok))
cat(sprintf("  Per-basin majority: PASS=%d, FAIL=%d, NA=%d\n", n_pass, n_fail, n_na))

# List FAIL diseases
fail_diseases <- atlas[direction_ok == FALSE, disease_family]
cat(sprintf("\n  FAIL diseases (%d): %s\n", length(fail_diseases), paste(fail_diseases, collapse = ", ")))

# ============================================================
# 4. Regenerate atlas heatmap
# ============================================================
cat("\n=== Regenerating atlas heatmap ===\n")

sm <- fread("summary/atlas_summary.tsv")
sm[, direction_label := ifelse(
  is.na(direction_ok), "N/A",
  ifelse(direction_ok == TRUE, "PASS", "FAIL")
)]
sm[, hr_label := sprintf("%.3f", HR)]
sm[, disease_family := factor(disease_family,
     levels = rev(sort(unique(disease_family))))]

g <- ggplot(sm, aes(x = outcome, y = disease_family,
                     fill = direction_label)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = hr_label), size = 3) +
  scale_fill_manual(
    values = c(PASS = "#4CAF50", FAIL = "#F44336", "N/A" = "gray80"),
    name = "Per-basin majority"
  ) +
  labs(title = "Energy Cox Sanity Check Atlas (Landmark Design)",
       subtitle = "Overall HR; color = per-basin majority criterion (HR<1 in ≥half basins)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))

ggsave("plot/atlas_heatmap.png", g, width = 8, height = 16, dpi = 200)
cat("  Written: plot/atlas_heatmap.png\n")

# Also save to SCRATCH for reference
ggsave(file.path(SCRATCH, "landmark_atlas_heatmap.png"), g, width = 8, height = 16, dpi = 200)

cat("\n=== Post-processing complete ===\n")
