## Identify diseases where energy is driven by disease-specific (non-generic) comorbidities
##
## 1. Define "generic" ICD-10 codes = those appearing in >50% of diseases' best PC top-10
## 2. For each disease: fraction of best PC loading weight from generic vs specific codes
## 3. Rank by specificity (low generic fraction = energy driven by non-standard codes)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "summary"

loadings <- fread(file.path(out_dir, "best_pc_loadings.tsv"))
best_pc <- fread(file.path(out_dir, "pc_energy_best.tsv"))

n_diseases <- length(unique(loadings$disease))

# --- Step 1: identify generic codes ---
code_freq <- loadings[, .N, by = icd10][order(-N)]
cat("=== ICD-10 code frequency in best PC top-10 loadings ===\n")
print(code_freq[1:20])

generic_threshold <- n_diseases * 0.5
generic_codes <- code_freq[N >= generic_threshold, icd10]
cat(sprintf("\nGeneric codes (>50%% of diseases = >%d): %d codes\n",
            as.integer(generic_threshold), length(generic_codes)))
cat(paste(generic_codes, collapse = ", "), "\n\n")

# --- Step 2: for each disease, compute generic vs specific loading weight ---
loadings[, abs_loading := abs(loading)]
loadings[, is_generic := icd10 %in% generic_codes]

disease_spec <- loadings[, .(
  total_weight = sum(abs_loading^2),
  generic_weight = sum(abs_loading[is_generic]^2),
  specific_weight = sum(abs_loading[!is_generic]^2),
  n_generic = sum(is_generic),
  n_specific = sum(!is_generic),
  top_specific_codes = paste(icd10[!is_generic][1:min(5, sum(!is_generic))],
                              collapse = ";")
), by = disease]

disease_spec[, generic_frac := generic_weight / total_weight]
disease_spec[, specific_frac := 1 - generic_frac]

# Merge with best PC info
disease_spec <- merge(disease_spec, best_pc[, .(disease, pc, cor_energy,
                                                 mort_pass)],
                      by = "disease")
disease_spec <- disease_spec[order(generic_frac)]

fwrite(disease_spec, file.path(out_dir, "pc_specificity.tsv"), sep = "\t")

cat("=== Disease ranking by specificity (low generic_frac = more specific) ===\n")
print(disease_spec[, .(disease, pc, generic_frac, specific_frac,
                        n_specific, mort_pass, top_specific_codes)], nrows = 60)

# --- Step 3: show top-10 specific loadings for exceptional diseases ---
exceptional <- disease_spec[generic_frac < 0.5, disease]
cat(sprintf("\n=== Exceptional diseases (generic_frac < 0.5): %d ===\n",
            length(exceptional)))

for (d in exceptional) {
  cat(sprintf("\n--- %s (best PC: %s, generic_frac: %.2f) ---\n",
              d, disease_spec[disease == d, pc],
              disease_spec[disease == d, generic_frac]))
  spec_load <- loadings[disease == d & !is_generic][order(-abs_loading)]
  if (nrow(spec_load) > 0) {
    print(spec_load[, .(rank, icd10, loading)])
  }
}

# --- Plot ---
disease_spec[, disease_f := factor(disease, levels = rev(disease))]
disease_spec[, pass_label := ifelse(mort_pass, "PASS", "FAIL")]

g <- ggplot(disease_spec, aes(x = generic_frac, y = disease_f,
                               fill = pass_label)) +
  geom_col() +
  geom_vline(xintercept = 0.5, linetype = 2, color = "gray30") +
  scale_fill_manual(values = c(PASS = "#1565C0", FAIL = "#C62828")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Generic vs specific: what fraction of energy is driven by age-related codes?",
       subtitle = "Low = energy driven by disease-specific comorbidities",
       x = "Generic code fraction (loading weight)",
       y = NULL, fill = "Mortality") +
  theme_minimal(base_size = 10)
ggsave(file.path(out_dir, "pc_specificity_ranking.png"),
       g, width = 11, height = 14, dpi = 150)

cat("\nDone.\n")
