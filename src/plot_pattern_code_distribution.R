## Distribution of ICD-10 codes per 128 pattern at multiple frequency thresholds
##
## For each disease, reads state_code_enrichment.tsv and counts how many codes
## each pattern has at thresholds 5%, 10%, 20%, 50%.
##
## Usage:
##   Rscript src/plot_pattern_code_distribution.R <disease> <out_dir>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
disease <- args[1]
out_dir <- args[2]

cat(sprintf("=== Pattern code distribution: %s ===\n", disease))

dt <- fread(file.path("plot", disease, "state_code_enrichment.tsv"))

thresholds <- c(0.05, 0.10, 0.20, 0.50)
threshold_labels <- c("5%", "10%", "20%", "50%")

results <- rbindlist(lapply(seq_along(thresholds), function(i) {
  th <- thresholds[i]
  per_pat <- dt[freq_in >= th, .(n_codes = uniqueN(code)),
                by = .(pattern_id, state_label)]
  # Patterns with 0 codes at this threshold
  all_pats <- unique(dt[, .(pattern_id, state_label)])
  per_pat <- merge(all_pats, per_pat, by = c("pattern_id", "state_label"), all.x = TRUE)
  per_pat[is.na(n_codes), n_codes := 0L]
  per_pat[, threshold := threshold_labels[i]]
  per_pat
}))

results[, threshold := factor(threshold, levels = threshold_labels)]

g <- ggplot(results, aes(x = n_codes)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~threshold, scales = "free", labeller = labeller(
    threshold = function(x) paste0("Freq >= ", x, " of patients"))) +
  labs(title = sprintf("%s: ICD-10 codes per pattern (128 patterns) at different thresholds",
                        disease),
       subtitle = "How many disease codes characterize each Ising pattern?",
       x = "Number of ICD-10 codes per pattern",
       y = "Number of patterns") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

png_path <- file.path(out_dir, "pattern_code_distribution.png")
ggsave(png_path, g, width = 12, height = 8, dpi = 200, bg = "white")
cat(sprintf("Written: %s\n", png_path))

# Summary stats
cat("\n=== Summary ===\n")
for (th in threshold_labels) {
  sub <- results[threshold == th]
  cat(sprintf("  %s: median=%.0f, mean=%.1f, min=%d, max=%d\n",
              th, median(sub$n_codes), mean(sub$n_codes),
              min(sub$n_codes), max(sub$n_codes)))
}

cat("Done.\n")
