## Regenerate FigS3: within-disease basin-to-basin Jaccard similarity
##
## Reads verified jaccard_all_basins.tsv (56 within-disease pairs).
## Basin count label = n_distinct(basin_id > 0), NOT max(basin_id).
## Reference line = mean of all between-disease basin-level pairs.
##
## Output: plot/Figures/Supplementary/FigS3_Fig04_within_disease_basin_similarity.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base <- "/home/koki/dev/DeSC-DiseaseLandscapeAtlas"
setwd(base)

jac <- fread("summary/jaccard_all_basins.tsv")
cat(sprintf("Total pairs: %d\n", nrow(jac)))

# Within-disease pairs
within <- jac[disease1 == disease2]
cat(sprintf("Within-disease pairs: %d\n", nrow(within)))

# Between-disease pairs for reference line
between <- jac[disease1 != disease2]
between_mean <- mean(between$jaccard)
cat(sprintf("Between-disease mean: %.4f (n=%d)\n", between_mean, nrow(between)))

# Count distinct basins per disease from the within-disease pairs
basin_counts <- list()
for (d in unique(within$disease1)) {
  ids <- unique(c(within[disease1 == d, basin1], within[disease1 == d, basin2]))
  basin_counts[[d]] <- length(ids)
}
bc_dt <- data.table(disease = names(basin_counts),
                    n_basins_distinct = unlist(basin_counts))

# Human-readable names
name_map <- c(
  acute_renal_failure = "Acute Renal Failure",
  anemia = "Anemia",
  anxiety = "Anxiety",
  aspiration_pneumonia = "Asp. Pneumonia",
  asthma = "Asthma",
  atrial_fibrillation = "Atrial Fibrillation",
  benign_prostatic_hyperplasia = "BPH",
  ckd = "CKD",
  delirium = "Delirium",
  dementia = "Dementia",
  depression = "Depression",
  diabetes = "Diabetes",
  dyslipidemia = "Dyslipidemia",
  epilepsy = "Epilepsy",
  heart_failure = "Heart Failure",
  hematologic_malignancy = "Hemat. Malignancy",
  hypertension = "Hypertension",
  inflammatory_bowel_disease = "IBD",
  liver_cancer = "Liver Cancer",
  malignancy = "Malignancy",
  malnutrition = "Malnutrition",
  osteoarthritis = "Osteoarthritis",
  osteoporosis = "Osteoporosis",
  parkinson = "Parkinson",
  peptic_ulcer = "Peptic Ulcer",
  pneumonia = "Pneumonia",
  prostate_cancer = "Prostate Cancer",
  rheumatoid_arthritis = "Rheumatoid Arthritis",
  schizophrenia = "Schizophrenia",
  sepsis = "Sepsis",
  sleep_disorders = "Sleep Disorders",
  thyroid_disorder = "Thyroid",
  urinary_tract_infection = "UTI",
  viral_hepatitis = "Viral Hepatitis"
)

within <- merge(within, bc_dt, by.x = "disease1", by.y = "disease")
within[, display_name := ifelse(disease1 %in% names(name_map),
                                name_map[disease1], disease1)]
within[, y_label := sprintf("%s (%d basins)", display_name, n_basins_distinct)]
within[, above := jaccard >= between_mean]

# Order diseases by mean within-disease Jaccard (descending, matching original layout)
disease_order <- within[, .(mean_jac = mean(jaccard)), by = y_label][order(mean_jac)]
within[, y_label := factor(y_label, levels = disease_order$y_label)]

# Verify pair counts
cat("\nPer-disease verification:\n")
for (d in unique(within$disease1)) {
  sub <- within[disease1 == d]
  n_b <- bc_dt[disease == d, n_basins_distinct]
  expected <- n_b * (n_b - 1) / 2
  ok <- nrow(sub) == expected
  cat(sprintf("  %s: %d basins, %d pairs (expected %d) %s\n",
              d, n_b, nrow(sub), expected, ifelse(ok, "OK", "MISMATCH")))
}

g <- ggplot(within, aes(x = jaccard, y = y_label, color = above)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_vline(xintercept = between_mean, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  annotate("text", x = between_mean - 0.005, y = nrow(disease_order) - 0.5,
           label = sprintf("Between-disease\nmean = %.3f", between_mean),
           hjust = 1, vjust = 1, size = 3.5, fontface = "italic") +
  scale_color_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C"),
                     labels = c("TRUE" = "Above between-disease mean",
                                "FALSE" = "Below between-disease mean"),
                     name = NULL,
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(limits = c(0.3, 1.0)) +
  labs(title = "Within-disease Basin Similarity",
       x = "Jaccard Similarity (within-disease basin pairs)",
       y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    legend.position.inside = c(0.70, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray90")
  )

out_path <- "plot/Figures/Supplementary/FigS3_Fig04_within_disease_basin_similarity.png"
ggsave(out_path, g, width = 10, height = 14, dpi = 200, bg = "white")
out_pdf <- sub("\\.png$", ".pdf", out_path)
ggsave(out_pdf, g, width = 10, height = 14, bg = "white")
cat(sprintf("\nSaved: %s\n", out_path))
cat(sprintf("Saved: %s\n", out_pdf))

cat(sprintf("\nTotal dots plotted: %d (expected 56)\n", nrow(within)))
